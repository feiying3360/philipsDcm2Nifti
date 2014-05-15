#!/usr/bin/env python
'''
A Philips MR multiframe dicom to nifti converter

MAIN LIMITATIONS:
Uses main PC memory to reformat data (could be an issue with large data files, eg fmri)

Created on 22 Jan 2013
@date: $Date: 2013-12-19 14:10:06 +0000 (Thu, 19 Dec 2013) $
@author: gfagiolo

The MIT License (MIT)

Copyright (c) 2013 Gianlorenzo Fagiolo

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


'''

#Standard Python Library
import os
import sys
from glob import glob
from optparse import OptionParser
from string import Template
import logging
from philips_dcm import PhilipsMultiFrameDcm
from philips_dcm.helpers import is_mrspectroscopystorage
logging.basicConfig(level=logging.WARNING)
import gzip

#Libraries from the internet
import numpy

#Own written
from pynii import Nifti1Data  # @UnresolvedImport
import philips_dcm

__version__ = "0.6 ($Revision: 32 $)"

def modify_revision_tag(ast):
    return ast.replace('$','').replace(' ','').replace(':','').replace('Revision','.').replace('(','').replace(')','')

DEFAULT_DESCRIPTION = "dcm_extractor " + \
    modify_revision_tag(__version__) + \
    ' philips_dcm ' + modify_revision_tag(philips_dcm.__version__)

DEFAULT_TIME_DIM = 1000. #msecs

def set_nifti_description(img):
    img.setDescription(DEFAULT_DESCRIPTION + ' ' + modify_revision_tag(img.getDescription()))

def save_mrs_voxel_nifti(pd, spect_vox_name):
    """creates a single voxel nifti
    @type pd: PhilipsMultiFrameDcm
    @type spect_vox_name: str 
    """
    nii = Nifti1Data()
    nii.setAffine(pd.get_minimal_header().get_nifti_affine())
    nii.setData(100 * numpy.ones((1, 1, 1), dtype=numpy.int16))
    set_nifti_description(nii) 
    Nifti1Data.save(nii, spect_vox_name)
    

def dcm_2_nifti(fname, outname=None, 
                verbose=False, debug_mode=False, 
                dest_dir=None, add_sequence_description=False,
                save_metadata_file=True, only_sv_mrs=False):
    #set export filename
    if outname is None:
        outname = str(fname)
    elif not dest_dir is None:
        outname = os.path.join(dest_dir, os.path.split(outname)[1])
    

    if verbose:
        numpy.set_printoptions(precision=3, suppress=True)
        hbar = 80*'#'
        print hbar
        print fname
    if only_sv_mrs and not is_mrspectroscopystorage(fname):
        raise ValueError('Not an MRS acquisition (mrs_only option is on)')
    pd = PhilipsMultiFrameDcm(fname)
    if add_sequence_description:
        try:
            outname += '.' + pd.get_sequence_custom_description_txt()
        except Exception:
            logging.warning("couldn't determine series number")
    #save metadata
    if save_metadata_file:
        save_metadata(pd, outname)
    if pd.is_mrs():
        #this is an MRSpectroscopy acquisition
        #create a single voxel nifti
        if not outname.endswith('.nii.gz'):
            spect_vox_name = outname + '.nii.gz'
        else:
            spect_vox_name = outname
        save_mrs_voxel_nifti(pd, spect_vox_name)
        print 'SUCCESS: Saved %s' % (spect_vox_name)    
    else:
        img_data = pd.get_pixel_array() 
        hd = pd.get_minimal_header()
        if debug_mode:
            #save list of indexes with header
            ofile = open(fname + '.slice_indices.csv', 'wb')
            ofile.write(hd.get_indices_list().toCSV())
            ofile.close()
        order = hd.default_order()
        orslcs = hd.sort_slices(order)
        stackid = hd.get_index_stack_id()
        imtid = hd.get_index_image_type_mr()
        echid = hd.get_index_effective_echo_time()
        nostacks = len(hd.get_stack_orientation().keys())
        if hd.is_dti():
            #remove slices that belong to reconstructed/fitted ADC image
            excluded = map(lambda x:x.get_frame_no(), hd.get_dti_info_list().get_processed_data())
            if excluded:
                orslcs = orslcs.get_slices_excluding_frame_list(excluded)
                #warn user
                logging.warning("DTI acquisition, only saving acquired data (excluded processed data such as AD images)")
                
        if verbose:
            template = "{0:35}# {1:6}"
            print hbar
            print template.format('Index Name','Index Count')
            print hbar
            for x in zip(
                    pd.get_minimal_header().get_dimension_indices_labels(),
                    pd.get_minimal_header().get_dimension_indices_counts()):
                print template.format(*map(str,x))
            print hbar
            print template.format('Information','Value')
            print hbar
            print template.format('Total Frames', str(len(orslcs)))
            print template.format('Order', ','.join([hd.get_dimension_indices_labels()[i] for i in order]))
            print hbar
            
        #Split data if different image type MR occur (since there are different intensity scaling slope and intercept)
        #Subsequently split data into stacks (since there are different orientation/position information [one affine per stack])
        
        #loop over image type MR (i.e. magnitude, real, imaginary, phase)
        for imt in hd.get_intensity_scaling().keys():
            #Each image type image can have different rescale information
            #get intensity rescale info
            resc = hd.get_intensity_scaling()[imt]        
            if imt is None:
                #only one image type => get all ordered slices
                scSet = philips_dcm.FrameIndexList(orslcs)
            else:
                #select slices from this type only 
                scSet = orslcs.select_slices(imtid, imt)
            
            if verbose:
                #print out rescale information            
                print template.format('Image Type MR', str(imt))
                print template.format('Slope', str(resc.get_slope()))
                print template.format('Intercept', str(resc.get_intercept()))
                print hbar
                
            #loop over stacks       
            for st in hd.get_stack_orientation().keys():
                #each stack can have different FOV setting
                #get affine transformation for this stack
                aff = hd.get_nifti_affine_stack(st)
                if verbose:
                    #print affine tranformation for DICOM and NIFTI
                    print 'Stack %d'%(st)
                    print hbar
                    print 'NIFTI AFFINE:'
                    print aff[:3,:]
                    print 'DICOM AFFINE:'
                    print hd.get_affine_stack(st)[:3, :]
                    print hbar
                    
                if debug_mode:
                    #when debug mode is on, no files are actually saved
                    #skip next instructions
                    continue
                #only select slices from current stack [the last element of scSet is the frame number]
                selSet = scSet.select_slices(stackid, st)
                
                #base output filename
                if not outname.endswith('.nii.gz'):
                    bname = outname + '.nii.gz'
                else:
                    bname = str(outname)
                    
                #if more than 1 stack, add stack number to file name
                if nostacks > 1:
                    bname = bname.replace('.nii', '.s%d.nii'%(st))
                #if more than 1 image type, add image type number
                if not imt is None:
                    bname = bname.replace('.nii', '.it%d.nii'%(imt))
                
                #split data onto different files if different echo times are present
                if not echid is None:                
                    #loop over different echoes
                    for echono in hd.get_echo_numbers():
                        #separate echoes into different nifti files
                        #add to filename the echo number
                        nname = bname.replace('.nii', '.e%02d.nii'%echono)
                        reformat_and_save(nname, hd, img_data, selSet.select_slices(echid, echono), aff, resc)
                        print 'SUCCESS: Saved %s (stack %d, echo %d)'%(nname, st, echono)
                else:
                    reformat_and_save(bname, hd, img_data, selSet, aff, resc)
                    print 'SUCCESS: Saved %s (stack %d)'%(bname, st)


def save_metadata(pd_obj, outname):
    #save metadata
    metadatafile = gzip.GzipFile(outname.replace('.nii','').replace('.gz','') + '.metadata.gz', 'wb')
    metadatafile.write(philips_dcm.METADATA_DESCRIPTION)
    metadatafile.write(philips_dcm.dicomobj_to_str(pd_obj.get_dicom_obj(), only_non_private=False))
    metadatafile.close()

def reformat_and_save(bname, hd, img_data, selected_frames, aff, resc):
    """
    reformats image data, saves it to nifti 
    it also saves extra files if it is a dti acquisition (bvals, bvecs and a bash script to run FSL dtifit)
    """
    nslices = hd.get_number_of_slices()
    img = Nifti1Data()
    if hd.format_as_4D():
        #TODO: encode 4th dimension pixdim
        #Save as 4D
        rem = int(len(selected_frames) /  nslices)
        nshape =   (hd.get_no_rows(), hd.get_no_cols(), nslices, rem)                      
        data = img_data[selected_frames.get_frames_indices_only(), :, :].transpose((2, 1, 0)).reshape(nshape)
    else:
        data = img_data[selected_frames.get_frames_indices_only(), :, :].transpose((2, 1, 0))
    #add description information
    set_nifti_description(img)
    #set affine
    img.setAffine(aff)
    #set data
#     img.setData(numpy.int16(data))
    img.setData(data)
    #update slope and rescale
    img.setScaleSlope(resc.get_slope())     
    img.setScaleIntercept(resc.get_intercept())
    #set spatial units       
    img.setXYZunit('mm')
    if hd.format_as_4D() and not hd.is_dti():
        #Write pixdim and unit for temporal dimension
        #ASSUMPTIONS: time unit is repetition time
        t_dim = hd.get_repetition_time()
        if t_dim is None:
            #drop back to default value and warn user
            t_dim = DEFAULT_TIME_DIM
            logging.warning("Couldn't find Repetition Time using default value for time unit (1000 msecs)")
        img.setTemporalDimension(t_dim)
        img.setTunit('msec')
    #save nifti
    Nifti1Data.save(img, bname)    
    #save extra information if it is a DTI acquisition
    if hd.is_dti():
        #dti acquisition
        #save bvecs and bvals in FSL format
        dti_idx = (hd.get_index_diffusion_bvalue(), hd.get_index_diffusion_grad_orient())
        try:
            hd.get_dti_info_list().save_FSL_bval_bvecs(selected_frames, dti_idx, bname + '.bvals', bname + '.bvecs')
            #create a linux script to run FSL dtifit
            save_FSL_script(bname)
        except TypeError as e:
            logging.error( 'WARNING: Could not save dti extra files bvals, bvecs etc. (possibly because this is a derived data set)\n' + str(e))      
    return bname
    
def save_FSL_script(niiname):
    script_t = Template("""#!/bin/sh

#check if FSL is properly installed

if [ -z "$$FSLDIR" ]
then
   echo FSL is not properly installed/configured. Please set the FSLDIR environment variable
    exit 1
fi

echo processing file: $infile

#do brain extraction

bet $infile $brainfile -m -R
echo brain extraction completed: $brainfile, $brainmaskfile

#do eddy current correction

eddy_correct $infile $eddycorrectfile 0
echo eddy current correction completed: $eddycorrectfile

#do tensor model fit

dtifit --data=$eddycorrectfile --out=$outname --mask=$brainmaskfile --bvecs=$bvecs --bvals=$bvals
echo dti fit completed: ${outname}_\(FA, MD, Ln, Vn\)

#view results with fslview

fslview ${outname}_FA ${outname}_V1
""")
    script_name = niiname + '.fsl_dtifit.sh'
    fname = os.path.split(niiname)[1]
    bname = fname.replace('.nii','.brain.nii')
    subs = {'infile':fname, 'brainfile':bname,
            'brainmaskfile':bname.replace('.nii','_mask.nii'),
            'eddycorrectfile':fname.replace('.nii','.ec.nii'),
            'outname':fname + '.dtifit',
            'bvals':fname + '.bvals',
            'bvecs':fname + '.bvecs',            
            }
    ofile = open(script_name, 'wb')
    ofile.write(script_t.substitute(subs))
    ofile.close()

def inspect_path(fname, options, cdepth=0):
    if not options.recursive and cdepth>0:
        return None
    if os.path.isfile(fname):
        try:
            dcm_2_nifti(fname, verbose=options.verbose, 
                        debug_mode=options.debug, 
                        dest_dir=options.output_directory,
                        add_sequence_description=options.add_sequence_description,
                        save_metadata_file=options.save_metadata_file, only_sv_mrs=options.mrs_only)
        except ValueError as e:
            logging.error('%s skipping (%s)'%(fname, str(e)))
        except AttributeError as e:
            logging.error('%s skipping (%s)'%(fname, str(e)))    
    else:
        if fname.find('*')>-1:
            #if there's a wild card, remove it
            fname = os.path.join(os.path.split(fname)[:-1])
        flist = glob(os.path.join(fname, '*'))
        if len(flist) < 1:
            logging.error('No files found by ' + fname)
            return None
        for t in flist:
            if os.path.isfile(t):
                try:
                    dcm_2_nifti(t, verbose=options.verbose, 
                                debug_mode=options.debug, 
                                dest_dir=options.output_directory,                                
                                add_sequence_description=options.add_sequence_description,
                                save_metadata_file=options.save_metadata_file, only_sv_mrs=options.mrs_only)
                except ValueError as e:
                    logging.error('%s skipping (%s)'%(t, str(e)))
                except AttributeError as e:
                    logging.error('%s skipping (%s)'%(t, str(e)))
            else:
                inspect_path(t, options, cdepth + 1)
    return None
    

def multi_file_cli():
    """Input:
    args and options are from an OptionParse
Output:
    Saves nifti files    
    """
    parser = OptionParser(usage="""%prog [options] dicom_file (dicom_file ... )
Converts the dicom_file specified to nifti files""", version=__version__)
    parser.add_option("-v", dest='verbose', action="store_true", default=False, help="prints useful information")
    parser.add_option("-d", dest='debug', action="store_true", default=False, help="does not save files and saves csv files containing slice order information.")
    parser.add_option("-o", dest='output_directory', default=None, help="where to store the converted nifti files")
    parser.add_option("-r", dest='recursive', action="store_true", default=False, help="recursively inspect path [not default]")    
    parser.add_option("-m", dest='multiple_files', action="store_true", default=False, help="specify multiple dicom files on the command line [not default]")
    parser.add_option("-s", dest='add_sequence_description', action="store_true", default=False, help="adds a sequence description txt to the output filename [not default]")
    parser.add_option("-e", dest='save_metadata_file', action="store_true", default=False, help="saves a file with all the dicom metadata [not default]")
    parser.add_option("-z", dest='mrs_only', action="store_true", default=False, help="Only generates nifti for Single Voxel MRS")    
    (options, args) = parser.parse_args()
    if len(args) > 0:
        for fname in args:
            inspect_path(fname, options)
#            if fname.find('*')>-1:
#                flist = glob(fname)
#                if len(flist) < 1:
#                    print 'No files found by',fname
#                    continue
#                for t in flist:
#                    try:
#                        dcm_2_nifti(t, verbose=options.verbose, debug_mode=options.debug, dest_dir=options.output_directory)
#                    except ValueError as e:
#                        print 'ERROR: %s skipping (%s)'%(t, str(e))
#                    except AttributeError as e:
#                        print 'ERROR: %s skipping (%s)'%(t, str(e))
#            else:
#                try:
#                    dcm_2_nifti(fname, verbose=options.verbose, debug_mode=options.debug, dest_dir=options.output_directory)
#                except ValueError as e:
#                        print 'ERROR: %s skipping (%s)'%(fname, str(e))
#                except AttributeError as e:
#                        print 'ERROR: %s skipping (%s)'%(fname, str(e))
    else:
        parser.print_help()

def simple_cli_api():
    parser = OptionParser(usage="""%prog [options] dicom_file outfile
Converts the dicom_file to nifti files with name specifide by outfile""", version=__version__)
    parser.add_option("-v", dest='verbose', action="store_true", default=False, help="prints useful information")
    parser.add_option("-d", dest='debug', action="store_true", default=False, help="does not save files and saves csv files containing slice order information.")
    parser.add_option("-m", dest='multiple_files', action="store_true", default=False, help="specify multiple dicom files on the command line [not default]")
    parser.add_option("-e", dest='save_metadata_file', action="store_false", default=True, help="Does not save a file with all the dicom metadata [it saves it]")    
    (options, args) = parser.parse_args()
    if len(args) == 2:
        try:
            dcm_2_nifti(args[0], args[1], verbose=options.verbose, debug_mode=options.debug,
                        save_metadata_file=options.save_metadata_file)
        except ValueError as e:
                logging.error('%s skipping (%s)'%(args[0], str(e)))
        except AttributeError as e:
                logging.error('%s skipping (%s)'%(args[0], str(e)))
    else:
        parser.print_help()
    
if __name__ == '__main__':    
    #decide whether to use non-default multi file CLI or single file CLI
    if '-m' in sys.argv:
        multi_file_cli()
    else:
        simple_cli_api()
    