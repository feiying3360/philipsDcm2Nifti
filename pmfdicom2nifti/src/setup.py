'''
Created on 15 May 2014

@author: glf12
'''
from distutils.core import setup

setup(name='dicom2nifti',
      version='0.6',
      description='A simple Philips MR Multiframe dicom to nifti converter',
      author='Gianlorenzo Fagiolo',
      url='https://github.com/gianlo/pynii',
      license='License :: OSI Approved :: MIT License',
      requires= ['numpy', 'dicom', 'philips_dcm'],
      py_modules=['dicom2nifti'])
