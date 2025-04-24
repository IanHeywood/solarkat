#!/usr/bin/env python


import re
import numpy
#import sys, os
#import subprocess
import numpy as np
from astropy.time import Time
from casacore.tables import table
from astropy import units as u
from MSUtils.msutils import addcol
#from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy.coordinates import solar_system_ephemeris, EarthLocation, AltAz
from astropy.coordinates import get_body #, get_body_barycentric


def hms2deg(hms):
    '''
    Function to convert hms to degrees
    '''
    hms_angle = Angle(hms, unit='hour')
    return hms_angle.degree


def dms2deg(dms):
    '''
    Function to convert dms to degrees
    '''
    dms_angle = Angle(dms, unit='degree')
    return dms_angle.degree


def extract_scan_number(ms_scan):
    # Extract the scan number from the scan file name using re
    scan_number = re.search(r"scan_(\d+)\.ms", ms_scan).group(1)
    return int(scan_number)


def rename_model_data_column1(ms, oldname, newname):

    '''
    Rename the 'MODEL_DATA' column to 'MODEL_DATA_ORIGINAL' in a Measurement Set (MS) table.
    Parameters:
    ms (str): Path to the Measurement Set file.
    oldname (str): Name of the column to be renamed.
    newname (str): New name for the column.
    '''
    print("Processing {}...".format(ms))
    # open the MS table in read-write mode
    ms = table(ms, readonly=False)
    # rename the 'MODEL_DATA' column to 'MODEL_DATA_ORIGINAL'
    ms.renamecol(oldname, newname)
    #print(f"Column '{oldname}' renamed to '{newname}' in {ms}.")
    print("Column '{}' renamed to '{}' in {}.".format(oldname, newname, ms))

    # close the MS table
    ms.close()


def rename_model_data_column(ms, oldname, newname):
    """
    Rename the 'MODEL_DATA' column to 'MODEL_DATA_ORIGINAL' in a Measurement Set (MS) table.
    If 'MODEL_DATA_ORIGINAL' already exists, it will be removed and replaced.
    Parameters:
    ms (str): Path to the Measurement Set file.
    oldname (str): Name of the column to be renamed.
    newname (str): New name for the column.
    """
    print(f"Processing: {ms}...")

    try:
        # Open the MS table in read-write mode
        with table(ms, readonly=False) as ms:
            # Check if the new column already exists
            if newname in ms.colnames():
                print(f"Column '{newname}' already exists. Removing it before renaming.")
                ms.removecols(newname)  # Remove the existing column
            # Check if the old column exists
            if oldname not in ms.colnames():
                print(f"Column '{oldname}' does not exist. Rename aborted.")
                return
            # Rename the column
            ms.renamecol(oldname, newname)
            print(f"Successfully renamed column '{oldname}' to '{newname}' in {ms}.")
    except Exception as e:
        print(f"An error occurred while processing {ms}: {e}")


def rename_columns(ms_list, oldname, newname):
    '''
    Rename columns in a list of MS tables.
    '''
    for ms_name in ms_list:

        print("Renaming column '{}' to '{}' in {}...".format(oldname, newname, ms_name))

        '''
        Open the MS table in read-write mode
        '''
        ms = table(ms_name, readonly=False)

         # Rename the 'oldname' column to 'newname'
        if 'oldname' in ms.colnames():
            ms.renamecol('oldname', 'newname')
        # close the MS table
        ms.close()
        print("Renaming completed for {}.".format(ms_name))
    print("Rename column completed successfully.")


def get_scans_info(ms, output_file):

    '''
    Extracts the coordinates of the Sun from a measurement set and writes them to a file.
    Also writes the scan number and the name of the eventual per-scan MS

    Parameters:
    ms (str): Path to the measurement set file.
    outfile (str): Path to the output file.

    Returns:
    None
    '''

    def format_coords(ra0,dec0):
        c = SkyCoord(ra0*u.deg,dec0*u.deg,frame='fk5')
        #hms = c.ra.to_string(u.hour, precision=1, pad=True)
        # Format Dec without decimal seconds
        #dms = c.dec.to_string(u.deg, precision=1, pad=True) 
        #return hms, dms
        hms = str(c.ra.to_string(u.hour, pad=True, precision=1))
        dms = str(c.dec.to_string(u.deg, pad=True, precision=1))
        #dms = str(c.dec)
        return hms,dms

    # MeerKAT
    obs_lat = -30.71323598930457
    obs_lon = 21.443001467965008
    loc = EarthLocation.from_geodetic(obs_lat, obs_lon)

    fieldtab = table(f'{ms}::FIELD')
    phase_dir = fieldtab.getcol("PHASE_DIR")
    fieldtab.close()

    # Extract the RA and Dec values from the PHASE_DIR column
    ra_orig = phase_dir[0,0,0]
    dec_orig = phase_dir[0,0,1]

    ra_orig_hms, dec_orig_dms = format_coords(ra_orig,dec_orig)

    maintab = table(ms)
    scans = list(numpy.unique(maintab.getcol('SCAN_NUMBER')))
    lines = []

    print(f"Determining per-scan solar coordinates from {ms}")
    for scan in scans:
        print(f"Processing scan {scan}")
        scan_ms = ms.replace('.ms',f'_scan_{scan}.ms')
        subtab = maintab.query(query='SCAN_NUMBER==' + str(scan))
        t_scan = numpy.mean(subtab.getcol('TIME'))
        t = Time(t_scan / 86400.0, format='mjd')
        subtab.close()

        with solar_system_ephemeris.set('builtin'):
            sun = get_body('Sun', t, loc)
            sun_ra = sun.ra.value
            sun_dec = sun.dec.value
            sun_ra_hms, sun_dec_dms = format_coords(sun_ra, sun_dec)
            lines.append(f"{scan_ms} {scan} {sun_ra_hms} {sun_dec_dms} {ra_orig_hms} {dec_orig_dms}")

    maintab.close()

    with open(output_file, 'wt') as f:
        for line in lines:
            f.write(line + '\n')

    print("Per scan information extracted and saved to {}.".format(output_file))
    print("Per-row format is scan_ms scan_number sun_ra_hms sun_dec_dms ra_orig_hms dec_orig_dms")


def read_scans_info(input_file):
    '''
    Open the scan info text file written by get_scans_info and return its contents
    '''
    scans_info = []
    with open(input_file, 'r') as file:
        for line in file:
            opms, scan, sun_ra, sun_dec, orig_ra, orig_dec = line.strip().split() # Assuming RA and Dec are separated by a space
            scans_info.append((opms, scan, sun_ra, sun_dec, orig_ra, orig_dec))
    return scans_info


def create_ds9_region_from_file(input_file, output_dir):
    """
    Function to create a DS9 region file for each coordinate in the scan info file
    """
    # #Open the MS table
    # tab = table(ms, readonly=True)
    # # Extract the scan numbers
    # scan_numbers = list(np.unique(tab.getcol('SCAN_NUMBER')))
    # # Close the MS table
    # tab.close()

    with open(input_file, 'r') as f:
        line = f.readline()
        while line:
            cols = line.split()
#            scan_ms = cols[0]
            scan_number = cols[1]
            ra = cols[2]
            dec = cols[3]
            print(f'Solar RA is {ra} and Dec is {dec} for {scan_number}')
            ra_deg = hms2deg(ra)
            dec_deg = dms2deg(dec)
            print(f'Creating DS9 region for scan {scan_number}')
            create_ds9_region(ra_deg, dec_deg, scan_number, output_dir)
            line = f.readline()
        f.close()
    print("DS9 region creation completed successfully.")


def create_ds9_region(ra, dec, scan_number, output_dir):
    # Function to create a DS9 region file
    sun_region = f"""# Region file format: DS9 CARTA 3.0.0-beta.3 
global color=green dashlist=8 3 width=3 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
circle({ra}, {dec}, 1188.0000") # color=green
"""
#The UHF band images of the sun has a radius of 1062.9118 and the L bands 1188.0000"
    with open(f"{output_dir}/sun-region-scan_{scan_number}.reg", 'w') as f:
        f.write(sun_region)


def add_column_to_ms(ms, colnames, likecol):
    """
    Add columns to a single measurement file
    """
    success = False
    try:
        print("Opening {}...".format(ms))
        tb = table(ms, readonly=False)
    except Exception as e:
        print("Error: {}".format(e))
        return success
    
    for colname in colnames:
        if colname not in tb.colnames():
            """
            Get column description from column 'like_col'
            """
            desc = tb.getcoldesc(likecol)
            desc[str('name')] = str(colname)
            desc[str('comment')] = str(desc['comment'].replace(" ", "_"))
            dminfo = tb.getdminfo(likecol)
            dminfo[str("NAME")] = "{}-{}".format(dminfo["NAME"], colname) 
            print("Adding column '{}' to {}...".format(colname, ms))       
            tb.addcols(desc, dminfo)
            success = True
    print("Columns added to {} successfully.".format(ms))
    tb.close()
    return success


def copy_model_data_to_model_data_sun(ms, ms_list, copycol, tocol):
    print(f"Copying model solar visibilities back to original MS")

    # Open the main MS table as maintab
    maintab = table(ms, readonly=False)

    # Get the unique scan numbers from the SCAN_NUMBER column in maintab
    scans = list(np.unique(maintab.getcol('SCAN_NUMBER')))

    for scan_number in scans:

        # Find the corresponding MS scan file for the current scan number
        ms_scan = next((ms_scan for ms_scan in ms_list if extract_scan_number(ms_scan) == scan_number), None)

        if ms_scan:
            # Open the MS scan as the source_subtab
            source_subtab = table(ms_scan, readonly=True)

            # Query the target_subtab to get the rows for the current scan number
            target_rows = maintab.query(query='SCAN_NUMBER==' + str(scan_number))

            # Get the MODEL_DATA from the source_subtab
            source_model_data = source_subtab.getcol(copycol)

            # Get the shape of the source_model_data
            source_shape = source_model_data.shape

            # Check if the target_subtab has the tocol column
            if tocol not in target_rows.colnames():
                # Add the tocol column to the target_subtab with the source_shape
                target_rows.addcols(
                    columns={tocol: {'datatype': 'complex', 'shape': source_shape}}
                )

            # Get the MODEL_DATA_SUN from the target_subtab
            target_model_data_sun = target_rows.getcol(tocol)

            # Update the target_model_data_sun with the source_model_data
            target_model_data_sun[:] = source_model_data

            # Put the updated MODEL_DATA_SUN back into the target_subtab
            target_rows.putcol(tocol, target_model_data_sun)

            # Flush the changes to disk
            maintab.flush()

            # Close the source_subtab
            source_subtab.close()
            print(f"   Copied {copycol} from {ms_scan} to {tocol} in {ms}")

    # Close the maintab
    maintab.close()

    print(f"Stored solar model visibilities in {tocol} in {ms}")


