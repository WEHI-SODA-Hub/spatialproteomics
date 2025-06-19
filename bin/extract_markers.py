#!/usr/bin/env python
'''
Module      : extract_markers
Description : Extracts markers from COMET OME-TIFF and writes comma-delimited
              output containing markers, exposure times and background channels
              compatible with input expected by mcmicro's
              background_subtraction tool.
Copyright   : (c) WEHI SODA Hub, 2025
License     : MIT
Maintainer  : Marek Cmero (@mcmero)
Portability : POSIX
'''

import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import typer
from typing import Annotated
from tifffile import TiffFile


def main(tiff_path: Annotated[Path, typer.Argument(
        help="Path to the TIFF input file.")]):
    """
    Extracts markers, exposure times, and background channels from COMET
    OME-TIFF
    """
    with TiffFile(tiff_path) as tiff:
        ome_metadata = tiff.ome_metadata

    root = ET.fromstring(ome_metadata)
    ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}

    # Extract channel name and index
    channels = root.findall('.//ome:Channel', ns)
    channel_meta = []
    for ch in channels:
        name = ch.attrib.get('Name')
        channel_index = ch.attrib.get('ID')
        channel_meta.append({'index': channel_index, 'marker_name': name})

    # Extract exposure times
    planes = root.findall('.//ome:Plane', ns)
    for pl in planes:
        c_index = int(pl.attrib.get('TheC'))
        exposure_time = pl.attrib.get('ExposureTime')
        # the planes are in the same order as the channels
        # but channel ID is not the same as this index
        channel_meta[c_index]['exposure'] = float(exposure_time)

    # Extract background channels
    channel_privs = root.findall('.//ome:ChannelPriv', ns)
    for idx, priv in enumerate(channel_privs):
        # just to be safe, let's make sure the ID is the same
        # (but it should be in the same order)
        channel_id = priv.attrib.get('ID')
        for cm in channel_meta:
            if cm['index'] == channel_id:
                background = priv.attrib.get('FluorescenceChannel')
                if background == cm['marker_name']:
                    # make it blank (background same as marker)
                    channel_meta[idx]['background'] = ""
                else:
                    channel_meta[idx]['background'] = background
                break

    # Format output and write to stdout
    df = pd.DataFrame(channel_meta)
    df = df[['marker_name', 'background', 'exposure']]
    df['exposure'] = df['exposure'].astype(float)
    df['remove'] = ''
    df.to_csv(sys.stdout, index=False)


if __name__ == "__main__":
    typer.run(main)
