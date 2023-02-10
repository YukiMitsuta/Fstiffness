#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright 2023/01/31 Yuki Mitsuta (mitsutay[at]omu.ac.jp)
#
# Distributed under terms of the MIT license.

import os
import numpy as np

def getFstiffness(Bandwidth,kappa,imagelist,Eholo):
    """
    function to calculate the stress of stiffness 
    BandWidth : the band width of stiffness term. 
                The larger this parameter is, the higher the stress of stiffness.
    kappa :     the spring constant of NEB (float or list of float)
                if the kappa is change in each connection (like as EW-NEB method) 
                plese prepaire the list of the spring constant.
    image :     the list of the coordinates of images
                the instance of each coordinates are numpy.array
    Eholo :     tangent vector of each images
    return
        Fstlist: the list of the stress of stiffness
    """
    Nimage = len(imagelist)
    Eholo_vertical, findvertlist = calcEholo_vert(imagelist,Eholo)
    if not any(findvertlist):
        return False
    innerimagelist = []
    outerimagelist = []
    for k in range(Nimage):
        innerimagelist.append(imagelist[k] - Eholo_vertical[k]*Bandwidth*0.5)
        outerimagelist.append(imagelist[k] + Eholo_vertical[k]*Bandwidth*0.5)
    Fstlist = [np.zeros(len(imagelist[k]))]
    for k in range(1, Nimage-1):
        if isinstance(kappa,list):
            Fst = Fstiffness_k(k, kappa[k], kappa[k+1], innerimagelist,outerimagelist, Eholo_vertical)
        else:
            Fst = Fstiffness_k(k, kappa, kappa, innerimagelist,outerimagelist, Eholo_vertical)
        Fstlist.append(Fst)
    return Fstlist
def Fstiffness_k(k,kappa0,kappa1,innerimagelist,outerimagelist,Eholo_vertical):
    """
    to calculate the stress of stiffness for each image
    """
    v0inner = innerimagelist[k]- innerimagelist[k+1]
    tdelta0inner = np.linalg.norm(v0inner)
    v0outer = outerimagelist[k]- outerimagelist[k+1]
    tdelta0outer = np.linalg.norm(v0outer)
    tdeltadelta = tdelta0inner - tdelta0outer
    Vvert = -Eholo_vertical[k]*tdeltadelta*kappa0*0.5

    v1inner = innerimagelist[k]- innerimagelist[k-1]
    tdelta1inner = np.linalg.norm(v1inner)
    v1outer = outerimagelist[k]- outerimagelist[k-1]
    tdelta1outer = np.linalg.norm(v1outer)
    tdeltadelta = tdelta1inner - tdelta1outer
    Vvert -= Eholo_vertical[k]*tdeltadelta*kappa1*0.5
    return Vvert
def calcEholo_vert(imagelist,Eholo):
    """
    function to calculate the verical vector of the tangent vector
    """
    calcEholovertTh = 0.001
    Nimage = len(imagelist)
    Eholo_vert = [None for _ in range(Nimage)] 
    findvertlist = [False for _ in range(Nimage)] 
    whileN = 0
    while whileN < 1000:
        whileN += 1
        if all(findvertlist):
            break
        for k in range(Nimage):
            if findvertlist[k]:
                continue
            tau = Eholo[k]
            if k == 0 and findvertlist[1]:
                if findvertlist[1]:
                    tau_nei = Eholo_vert[1]
                    a = -np.dot(tau_nei,tau)/np.dot(tau,tau)
                    Eholo_vert[k] = a*tau+tau_nei
                    findvertlist[k] = True
            elif k == Nimage-1:
                if findvertlist[-2]:
                    tau_nei = Eholo_vert[-2]
                    a = -np.dot(tau_nei,tau)/np.dot(tau,tau)
                    Eholo_vert[k] = a*tau+tau_nei
                    findvertlist[k] = True
            else:
                v1 = imagelist[k-1]-imagelist[k]
                v1 = v1/np.linalg.norm(v1)
                v1taudot = np.abs(np.dot(v1,tau))
                v2 = imagelist[k+1]-imagelist[k]
                v2 = v2/np.linalg.norm(v2)
                v2taudot = np.abs(np.dot(v2,tau))
                if 1.0-calcEholovertTh<v1taudot and 1.0-calcEholovertTh<v2taudot:
                    if findvertlist[k-1]:
                        tau_nei = Eholo_vert[k-1]
                        a = -np.dot(tau_nei,tau)/np.dot(tau,tau)
                        Eholo_vert[k] = a*tau+tau_nei
                        findvertlist[k] = True
                    elif findvertlist[k+1]:
                        tau_nei = Eholo_vert[k+1]
                        a = -np.dot(tau_nei,tau)/np.dot(tau,tau)
                        Eholo_vert[k] = a*tau+tau_nei
                        findvertlist[k] = True
                elif calcEholovertTh <= v1taudot:
                    a = -np.dot(v2,tau)/np.dot(v1,tau)
                    Eholo_vert[k] = a*v1 + v2
                    Eholo_vert[k] /= np.linalg.norm(Eholo_vert[k])
                    findvertlist[k] = True
                elif calcEholovertTh <= v2taudot:
                    a = -np.dot(v1,tau)/np.dot(v2,tau)
                    Eholo_vert[k] = a*v2 + v1
                    Eholo_vert[k] /= np.linalg.norm(Eholo_vert[k])
                    findvertlist[k] = True
                elif v1taudot < calcEholovertTh and v2taudot < calcEholovertTh:
                    Eholo_vert[k] = v1
                    findvertlist[k] = True
                else:
                    print("Error in calcEholovertTh")
                    exit()
    for k in range(1,Nimage):
        vbefore = Eholo_vert[k-1]
        v = Eholo_vert[k]
        if np.dot(vbefore,v) < 0.0:
            Eholo_vert[k] *= -1
    return Eholo_vert, findvertlist

