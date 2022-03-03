#!/usr/bin/env python3

# script to explore all folders where an edgelist is stored,
# compute the ML parameters and compute the projection

# # Dependencies
import csv      # for reading the edgelists
import os       # for os Ops
import NEMtropy # for automated analysis of the projected networks
import pickle   # for storage
import logging  # for logging

datapath = "/path/to/Twitter_IOps"

def myfun(folder, filename):
    """
        myfun(folder, filename)

    build model from edgelist, compute the ML parameters and store the pickled model
    """
    # load up the data
    logging.info("working on {}".format(os.path.join(folder, filename)))
    edges = edgeloader(os.path.join(folder, filename))
    
    # build and solve model
    model = NEMtropy.BipartiteGraph(edgelist=edges)
    model.solve_tool(method="fixed-point", initial_guess="degrees")
    # make projections
    model.compute_projection(threads_num=10)

    # store the results
    with open(os.path.join(folder, "model.pkl"), "wb") as outp:
        pickle.dump(model, outp, -1)
    logging.info("finished!".format(os.path.join(folder, filename)))

def edgeloader(p):
    """
    edgeloader

    load up an edgelist
    """
    res = []
    with open(p) as f:
        for line in f.readlines():
            s,d = line.split(sep=",")
            res.append((int(s), int(d)))

    return res

def list_files(dir):
    """
        list_files(dir)
    
    Go over all files and compute the pojected network.
    """
    for root, dirs, files in os.walk(dir):
        for direc in dirs:
            logging.info("working on folder: {}".format(direc))
            for rp, dd, df in os.walk(os.path.join(root,direc)):
                if "bipartite_edgelist.csv" in df:
                    # do something if edgelist exists
                    try:
                        myfun(rp, "bipartite_edgelist.csv")
                    except:
                        logging.warn("problem for {}".format(rp))
                else:
                    logging.warning("no edgelist found for {}".format(rp))

    logging.info("ALL DONE")

def main(path=datapath):
    logging.basicConfig(filename=os.path.join(path,'projector.log'), 
                        level=logging.INFO,
                        format='%(asctime)s - %(levelname)s:%(message)s')
    logging.info('Started')
    list_files(path)
    logging.info('Finished')


if __name__ == '__main__':
    main()
