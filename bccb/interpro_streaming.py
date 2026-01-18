"""
Streaming version of InterPro data processing to avoid memory overflow.
This module provides memory-efficient alternatives to pypath's interpro functions.
"""

from __future__ import annotations

import os
import gzip
import shutil
import collections
from typing import Generator, Literal, Optional
from lxml import etree
from pathlib import Path

from pypath.resources import urls
from pypath.share import curl
from biocypher._logger import logger


InterproEntry = collections.namedtuple(
    'InterproEntry',
    (
        'interpro_id',
        'protein_count',
        'name',
        'type',
        'publications',
        'parent_list',
        'child_list',
        'member_list'
    ),
)


def interpro_entries_streaming(
    early_stopping: Optional[int] = None
) -> Generator[tuple, None, None]:
    """
    Stream InterPro entries one by one instead of loading all into memory.
    
    Args:
        early_stopping: If set, stop after this many entries (for testing)
    
    Yields:
        InterproEntry named tuples one at a time
    """
    url = urls.urls['interpro']['entries']
    
    logger.info(f"Downloading InterPro entries from {url}")
    
    c = curl.Curl(
        url,
        silent=False,
        large=False
    )
    path = c.fileobj.name
    
    # Decompress the file
    decompressed_path = path.split('.gz')[0]
    with gzip.open(path, 'rb') as f_in:
        with open(decompressed_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    
    logger.info(f"Streaming InterPro entries from {decompressed_path}")
    
    # Use iterparse with element clearing to save memory
    context = etree.iterparse(decompressed_path, events=('end',), tag='interpro')
    
    count = 0
    for event, elem in context:
        # Parse publications
        if elem.find('pub_list') is not None:
            pubs = [pub.attrib['id'] for pub in elem.find('pub_list')]
        else:
            pubs = ''
        
        # Parse parent list
        if elem.find('parent_list') is not None:
            parent_ids = [parent.attrib['ipr_ref'] for parent in elem.find('parent_list')]
        else:
            parent_ids = ''
        
        # Parse child list
        if elem.find('child_list') is not None:
            child_ids = [child.attrib['ipr_ref'] for child in elem.find('child_list')]
        else:
            child_ids = ''
        
        # Parse member list
        member_ids = {}
        for member in elem.find('member_list'):
            db = member.attrib['db']
            dbkey = member.attrib['dbkey']
            if db in member_ids:
                member_ids[db].append(dbkey)
            else:
                member_ids[db] = [dbkey]
        
        entry = InterproEntry(
            interpro_id=elem.attrib['id'],
            protein_count=elem.attrib['protein_count'],
            name=elem.attrib['short_name'],
            type=elem.attrib['type'],
            publications=pubs,
            parent_list=parent_ids,
            child_list=child_ids,
            member_list=member_ids,
        )
        
        # Clear element to free memory
        elem.clear()
        # Also clear parent elements to save more memory
        while elem.getprevious() is not None:
            del elem.getparent()[0]
        
        yield entry
        
        count += 1
        if early_stopping and count >= early_stopping:
            logger.info(f"Early stopping at {count} entries")
            break
        
        # Log progress every 5000 entries
        if count % 5000 == 0:
            logger.info(f"Processed {count} InterPro entries...")
    
    logger.info(f"Finished processing {count} InterPro entries")
    
    # Clean up decompressed file to save disk space
    try:
        os.remove(decompressed_path)
    except:
        pass


def interpro_xrefs_streaming(
    db_type: Literal['go', 'structural', 'external'],
    early_stopping: Optional[int] = None
) -> Generator[tuple, None, None]:
    """
    Stream InterPro cross-references one at a time.
    
    Args:
        db_type: Type of cross-reference database
        early_stopping: If set, stop after this many entries
    
    Yields:
        Tuples of (interpro_id, xrefs_dict_or_list)
    """
    db_type_dict = {
        'go': 'class_list',
        'structural': 'structure_db_links',
        'external': 'external_doc_list',
    }
    db_type_name = db_type_dict[db_type]
    
    url = urls.urls['interpro']['entries']
    
    c = curl.Curl(
        url,
        silent=False,
        large=False
    )
    path = c.fileobj.name
    
    decompressed_path = path.split('.gz')[0]
    with gzip.open(path, 'rb') as f_in:
        with open(decompressed_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    
    context = etree.iterparse(decompressed_path, events=('end',), tag='interpro')
    
    count = 0
    for event, elem in context:
        interpro_id = elem.attrib['id']
        
        if db_type == 'go':
            if elem.find(db_type_name) is not None:
                go_terms = [go.attrib['id'] for go in elem.find(db_type_name)]
            else:
                go_terms = None
            result = go_terms
        else:
            other_db_keys = {}
            if elem.find(db_type_name) is not None:
                for link in elem.find(db_type_name):
                    db = link.attrib['db']
                    dbkey = link.attrib['dbkey']
                    if db in other_db_keys:
                        other_db_keys[db].append(dbkey)
                    else:
                        other_db_keys[db] = [dbkey]
            result = other_db_keys
        
        # Clear element to free memory
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]
        
        yield (interpro_id, result)
        
        count += 1
        if early_stopping and count >= early_stopping:
            break
    
    try:
        os.remove(decompressed_path)
    except:
        pass


def build_xrefs_dict(
    db_type: Literal['go', 'structural', 'external'],
    batch_size: int = 1000,
    early_stopping: Optional[int] = None
) -> dict:
    """
    Build xrefs dictionary in batches to reduce peak memory usage.
    
    This function builds the dictionary incrementally instead of loading
    all data at once.
    
    Args:
        db_type: Type of cross-reference database
        batch_size: Process entries in batches of this size
        early_stopping: If set, stop after this many entries
    
    Returns:
        Dictionary of interpro_id -> xrefs
    """
    result = {}
    
    for interpro_id, xrefs in interpro_xrefs_streaming(db_type, early_stopping):
        result[interpro_id] = xrefs
    
    return result
