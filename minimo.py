#!/usr/bin/env python3

from sys import exit
from os.path import isfile, join, split, splitext
import argparse
import netCDF4
import numpy as np
import logging

def pre_process_eq(eq):
    #this is for eq comming from lrs
    eq = np.unique(eq, axis = 0)
    eq = eq[eq[:, -1] != 0]
    # find the plane with biomass equal to zero
    non_zero_rows = np.any(eq[:, :-1], axis=1)
    eq = eq[non_zero_rows]
    return eq

def pre_process_aux_eq(eq):
    #this is for eq comming from lrs
    eq = np.unique(eq, axis = 0)
    eq = eq[eq[:, -1] != 0] 
    # find the plane with biomass equal to zero
    mask1 = np.any(eq[:, :-2], axis=1) | (np.all(eq[:, :-2] == 0, axis=1) & (eq[:, -1] != 0))
    eq = eq[mask1]

    # find the plane with aux equal to zero
    mask2 = np.any(eq[:, :-1], axis=1)
    eq = eq[mask2]
    return eq

def get_env(var,data):
    env_list = [None] * len(var)
    for idx, i in enumerate(var):
        compressed_data = data[i][:,:,:,:].compressed()
        env_list[idx] = compressed_data.reshape((compressed_data.shape[0], 1))
    return np.hstack(env_list)

def process_env_lrs(eq,env):

    eq_pos = eq[eq[:,-1]>0,:]
    eq_neg = eq[eq[:,-1]<0,:]

    # Initialize max_bio and min_bio
    max_bio = np.full((env.shape[0]), np.finfo(float).max)
    min_bio = np.full((env.shape[0]), np.finfo(float).min)

    if eq_pos.shape[0]>30:
        max_bio = np.full((env.shape[0]),np.finfo(float).max)
        for i in range(eq_pos.shape[0]):
            raw_bio = -((np.dot(eq_pos[i,1:-1],np.transpose(env)).transpose() \
                 + eq_pos[i,0])/eq_pos[i,-1]).transpose()
            max_bio = np.maximum(max_bio,raw_bio)
    else:
        if eq_pos.shape[0]>0:
            raw_bio = -((np.dot(eq_pos[:,1:-1],np.transpose(env)).transpose() \
                     + eq_pos[:,0])/eq_pos[:,-1]).transpose()
            max_bio = np.max(raw_bio, axis=0)
            print('Got max')
        else:
            max_bio = None

    if eq_neg.shape[0]>30:
        min_bio = np.full((env.shape[0]),np.finfo(float).max)
        for i in range(eq_neg.shape[0]):
            raw_bio = -((np.dot(eq_neg[i,1:-1],np.transpose(env)).transpose() \
                 + eq_neg[i,0])/eq_neg[i,-1]).transpose()
            min_bio = np.minimum(min_bio,raw_bio)
    else:
        raw_bio = -((np.dot(eq_neg[:,1:-1],np.transpose(env)).transpose() \
                 + eq_neg[:,0])/eq_neg[:,-1]).transpose()
        min_bio = np.min(raw_bio, axis=0)
    
    # test if max is < min
    # if max_bio != -1 and np.any(min_bio > max_bio):
    #     print('issue in computation')
    if max_bio is not None:
        print('verif this line')
        bio = np.maximum(max_bio,min_bio)
    else:
        bio = min_bio

    return bio

def get_plan(eq,env):
    # lrs is doing b A >= 0
    if eq.shape[0]>60:
        print("not implemented")
        return 0
    raw_bio = -((np.dot(eq[:,1:-1],np.transpose(env)).transpose() \
             + eq[:,0])/eq[:,-1]).transpose()
    min_bio = np.argmax(raw_bio, axis=0)
    return min_bio

def get_distance_to_stress(eq,env):
    # lrs is doing b A >= 0
    raw_bio = -((np.dot(eq[:,1:-1],np.transpose(env)).transpose() \
             + eq[:,0])/eq[:,-1]).transpose()
    min_bio = np.min(raw_bio, axis=0)
    bio_stress = raw_bio-min_bio
    return bio_stress

def read_lrs_hrep( filename , decimals = 8 ):
    '''
    This function should return the array that is in the file:
*lrs:lrslib v.6.2 2016.3.28(64bit,lrsgmp.h gmp v.6.1)
*Copyright (C) 1995,2016, David Avis   avis@cs.mcgill.ca 
*Input taken from file phaeo_photeuk_1034_comp_ready_lrs.ext
*volume
H-representation
begin
***** 4 rational
***
end
*Volume= 803834082439834112684634041415756860061592770497449/5846006549323611672814739330865132078623730171904 
*Totals: facets=10 bases=10
*Dictionary Cache: max size= 5 misses= 0/9   Tree Depth= 4
*lrs:lrslib v.6.2 2016.3.28(64bit,lrsgmp.h)
*0.002u 0.002s 3580Kb 2 flts 0 swaps 208 blks-in 8 blks-out
    '''
    with open(filename,'r') as f:
        lines = f.read().splitlines()
        first_line = 0
        is_first_line = (lines[first_line] == "H-representation")
        while not is_first_line:
            first_line+=1
            is_first_line = (lines[first_line] == "H-representation")
        last_line = len(lines) - 1
        is_last_line = (lines[last_line] == "end")
        while not is_last_line:
            last_line-=1
            is_last_line = (lines[last_line] == "end")
        first_used_line = 1
        is_first_line = (lines[first_line+first_used_line] == "begin")
        while not is_first_line:
            first_used_line+=1
            is_first_line = (lines[first_line+first_used_line] == "begin")
        used_lines = lines[first_line+first_used_line:last_line]
        # first line should have "begin"
        if used_lines[0] != "begin":
            print("warning not the expected format")
            print(used_lines[0])
        # next line should have "***** 4 rational"
        component_number = used_lines[1].split()[1]
        if not component_number.isdigit():
            print("Error in trying to retrieve the number of component")
            return -1
        component_number = int(component_number)
        print("We are looking at a file of dimension ",component_number)
        # next lines should have the numbers
        line_number = len(used_lines) - 2
        np_hrep = np.zeros((line_number,component_number))
        for i in range(2,len(used_lines)):
            numbers = used_lines[i].split()
            # print(numbers)
            for j in range(component_number):
                np_hrep[i-2,j] = float(round(Fraction(numbers[j]),decimals))
    return np_hrep

if __name__ == "__main__":
    '''
    TODO
    '''

    logging.basicConfig(level = logging.INFO)
    logger = logging.getLogger()
    parser = argparse.ArgumentParser()
    parser.add_argument('environmental_fluxes', 
        help='File with environmental fluxes in HDF4 format (.nc), temporal unit is second')
    parser.add_argument('niche_file', 
        help='File with the niche description in .ine (h-format), '
        'see lrslib for detail on the format. Unit of flux should be in molXX.gDW-1.h-1')
    parser.add_argument('flux_name', 
        help='File where flux names are stored, one per line, '
        'corresponding to the niche fluxes.')
    parser.add_argument('carbon', 
        help='Carbon composition of the considered system in mole.' 
        ' In a GSM, this is the stoichiometry of C in the biomass reaction.')

    parser.add_argument('-d','--directory', 
        help='Directory where result will be written, if not provided '
        'nothing will be written.')

    parser.add_argument('-s','--stress', 
        help='Nutrient on which to compute the stress')
    parser.add_argument('-a','--aux', 
        help='Auxiliary metabolite niche file computation')
    # parser.add_argument('-u','--unit', 
    #     help='Unit file that should have unit for XXX') # TODO
    # add a min of objective function
    # add a logger level

    args = parser.parse_args()

    # check for arguments

    filename_env = args.environmental_fluxes
    logger.info("Filename for the environment is " + str(filename_env) + ".")
    logger.info("Loading the environment..")

    # Read data
    data_env = netCDF4.Dataset(filename_env)

    # Look for every input in the dataset
    # TODO

    # Read niche file
    filename_niche = args.niche_file
    eq = read_lrs_hrep(filename_niche)

    # Read flux name file
    flux_file = args.flux_name
    with open(filename,'r') as f:
        fluxes = f.read().splitlines()
    for i in fluxes:
        try:
            data_env[i]
        except:
            logger.error("An exception occurred when retreiving fluxes " + i + 
                " in file " + filename_env)

    # check the size of the niche space with the number of line in the flux name file
    if eq.shape[0] > len(fluxes) + 3:
        logger.error("Niche dimension is " + str(eq.shape[0]) + 
            " while only " + str(len(fluxes)) + " fluxes are given. "
            "The fluxes in this file should be the same as the fluxes used "
            "for the niche computation.")
    if eq.shape[0] < len(fluxes) + 2:
        logger.error("Niche dimension is only " + str(eq.shape[0]) + 
            " while " + str(len(fluxes)) + " fluxes are given. "
            "The fluxes in this file should be the same as the fluxes used "
            "for the niche computation.")

    logger.info("Preprocessing..")
    p_eq = pre_process_eq(eq)

    env = get_env(fluxes,data_env)

    # Get carbon composition
    c_bio = args.carbon
    env *= 3600*c_bio

    logger.info("Computing..")
    bio = process_env_lrs(p_eq, env)

    # write bio
    if args.directory is not None:
        np.save(join(arg.directory,'growth.npy'),bio)

    # If asked compute stress
    if args.stress:
        bio_stress = get_distance_to_stress(p_eq,env)
        bio_stress = np.asarray(bio_stress,dtype=float)
        stress = bio_stress[fluxes.index(args.stress),:]
        # write stress
        if args.directory is not None:
            np.save(join(arg.directory,'stress.npy'),stress)

    # If asked compute auxiliary metabolite
    if args.aux is not None:
        # create env aux
        aux_env = np.hstack((env, bio.reshape((len(bio),1))))
        # read stress eq
        eq_aux_file = args.aux
        aux_eq = read_lrs_hrep(eq_aux_file)
        aux_p_eq = pre_process_aux_eq(aux_eq)
        # compute aux
        aux = process_env_lrs(aux_p_eq, aux_env)
        if args.directory is not None:
            np.save(join(arg.directory,'auxiliary.npy'),aux)

