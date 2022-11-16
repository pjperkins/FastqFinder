import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import datetime
import sys
import statistics

def process(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}

def readGz(fn = None):
    n = 4
    records = []
    with gzip.open(fn, 'r') as fh:
        lines = []
        for line in fh:
            lines.append(line.rstrip())
            if len(lines) == n:
                record = process(lines)
                #sys.stderr.write("Record: %s\n" % (str(record)))
                records.append(record)
                lines = []           
    return records
def read(fn = None):
    n = 4
    records = []
    with open(fn, 'r') as fh:
        lines = []
        for line in fh:
            lines.append(line.rstrip())
            if len(lines) == n:
                record = process(lines)
                #sys.stderr.write("Record: %s\n" % (str(record)))
                records.append(record)
                lines = []           
    return records

def fastqFinder(startpath = None):
    fils = []
    for path, currentDirectory, files in os.walk(startpath):
        for file in files:
            fils.append(os.path.join(path, file))
    fastqind = [i for i, s in enumerate(fils) if 'fastq' in s]
    fastqs = [fils[i] for i in fastqind]
    returndicts = []
    for fastq in fastqs:
        fil = fastq.split('/')
        fillen = len(fil) - 1
        fil = fil[fillen]
        ctime = os.path.getctime(fastq)
        ctime = datetime.datetime.fromtimestamp(ctime)
        size = os.path.getsize(fastq) / 1000
        if('.gz' in fil):
            data = readGz(fastq)
        if('.gz' not in fil):
            data = read(fastq)
        data = read(fastq)
        sequences = [d['sequence'] for d in data]
        lens = [len(d) for d in sequences]
        medlen = statistics.median(lens)
        returndict = {"path": fastq, 'file': fil, "size": size, "date": ctime, "length": medlen}
        returndicts.append(returndict)
    return returndicts


def plotLens(returndicts = None):
    lens = [d['length'] for d in returndicts]
    nams = [d['file'] for d in returndicts]
    df = pd.DataFrame(lens, nams)
    df.plot.bar(grid = True, color = '#607c8e')
    
    
def plotSizes(returndicts = None):
    sizes = [d['size'] for d in returndicts]
    nams = [d['file'] for d in returndicts]
    df = pd.DataFrame(sizes, nams)
    df.plot.bar(grid = True, color = '#607c8e')
    
    
def plotTimeline(returndicts = None):
    

