import time
import os
import argparse
import json
import lzma
import allium
import torch

import pickle

def main(args):
    """
    Basic call script for running cell migration inference using active brownian particle model
    """
    # Preparing simulations
    if args.simulate:
        if args.save_prob:
            print("saving data")
#         d = json.loads(args.thetakeys)        
        sim = allium.simulate.Sim(parameterFile=args.configfile,keys=args.thetakeys)

        obs = sim.simulate(args.theta)
        if args.save_prob == 1:
            picklefile = lzma.open(f"{args.outputfile}.pz", "wb")
            pickle.dump(obs, picklefile)
            
        ssvect, ssdata = allium.summstats.calculate_summary_statistics(
            obs,
            opts=args.summstatsopts,
            log=args.log,
            starttime=args.starttime,
            endtime=args.endtime,
        )
        if args.save_prob == 1:
            with lzma.open(f"{args.outputfile}_ss.pz", "wb") as f:
                pickle.dump(ssdata, f)
        return obs#, ssvect, ssdata

    if not os.path.exists(args.posteriorfolder):
        os.makedirs(args.posteriorfolder)
    if not os.path.exists(args.outputfolder):
        os.makedirs(args.outputfolder)
        os.makedirs(f"{args.outputfolder}/data")
#     d = json.loads(args.thetadict)
    thetamin = [float(t) for t in args.thetamin[: len(args.thetakeys)]]
    thetamax = [float(t) for t in args.thetamax[: len(args.thetakeys)]]
    prior = allium.utils.init_prior([thetamin, thetamax])
    posterior_opt = ["flow", "mix"]
    print("Preparing simulator")
    # Prepare simulation object
    sim = allium.simulate.Sim(
        keys=args.thetakeys,
        run=args.outputfile,
        parameterFile=args.configfile,
        ssopts=args.summstatsopts,
        framerate=args.framerate,
        log=args.log,
        test=args.test,
        folder=args.outputfolder,
        save_prob=args.save_prob,
        starttime=args.starttime,
        endtime=args.endtime,
        proposal=prior,
        num_simulations=args.nruns,
        batch_size=args.batch_size,
        num_workers=args.nprocs,
        useall=args.useall,
        maxN=args.maxN
    )

    # Running simulations
    print("Beginning simulation rounds:")
    tic = time.perf_counter()

    theta, x = sim.sample()

    # Save parameter/observable data
    picklefile = open(
        f'{args.outputfolder}data/{args.outputfile}_{time.strftime("%Y-%m-%d-%H-%M-%S", time.gmtime())}.p',
        "wb",
    )
    pickle.dump((theta, x), picklefile, )
    toc = time.perf_counter()
    print(f"Completed simulations in {toc - tic:0.4f} seconds")
    if len(args.starttime) > 1:
        if args.state == "confluent":
            x = x[:, :, 0].float()
        else:
            x = x[:, :, 1].float()

    toc = time.perf_counter()
    print(f"Completed simulating in {toc - tic:0.4f} seconds")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Running cell migration inference using active brownian particle models and neural network approximators"
    )
    # general options
    parser.add_argument("-c", "--configfile", default="include/config/simconfig.json")
    parser.add_argument("--sim", dest="simulate", action="store_true")
    parser.add_argument(
        "--test",
        dest="test",
        action="store_true",
        help="Testing summary statistics using previously saved file",
    )
    parser.add_argument("--no-test", dest="test", action="store_false")
    parser.add_argument(
        "--useall", dest="useall", action="store_true", help="Use all summary statistics data"
    )
    parser.add_argument(
        "--usepart", dest="useall", action="store_false", help="Use summarys of summary statistics data"
    )

    parser.add_argument(
        "--log", dest="log", action="store_true", help="Output_data to log_file"
    )
    parser.add_argument("--no-log", dest="log", action="store_false")
    parser.add_argument(
        "-s",
        "--save_prob",
        type=float,
        default=1,
        help="(float)\nProbability of saving file",
    )
    parser.add_argument(
        "-nruns",
        "--nruns",
        type=int,
        default=1,
        help="(int)\nNumber of simulations to run",
    )
    parser.add_argument(
        "-nprocs", "--nprocs", type=int, default=2, help="(int)\nNumber of cores to use"
    )
    parser.add_argument(
        "-batch", "--batch_size", type=int, default=1, help="(int)\nBatchsize"
    )
    # io data
    parser.add_argument(
        "-ofo",
        "--outputfolder",
        default="output/confluent/",
        help="(str)\nFolder for data",
    )
    parser.add_argument(
        "-ofi", "--outputfile", default="interaction", help="(str)\nFolder for data"
    )
    parser.add_argument(
        "-pfo",
        "--posteriorfolder",
        default="posteriors/",
        help="(str)\nFolder for posteriors",
    )
    parser.add_argument(
        "-pfi", "--pfile", default="posterior", help="(str)\nFilename for posteriors"
    )
    # simulation data (inc parameter and ranges)
    parser.add_argument(
        "-theta",
        "--theta",
        nargs="+",
        default=[],
        help="(list)\nList of parameter values to pass to simulation",
    )
    parser.add_argument(
        "-k",
        "--thetakeys",
        nargs="+",
        default=["factive", "pairatt", "tau", "phi"],
        help="(list)\nStrings for specific simulation parameters to update",
    )
    parser.add_argument(
        "-thetamin",
        "--thetamin",
        nargs="+",
        default=[0, 0, 10, 200],
        help="(list)\nList of lowerbound parameters values",
    )
    parser.add_argument(
        "-thetamax",
        "--thetamax",
        nargs="+",
        default=[0.2, 0.15, 400, 1200],
        help="(list)\nList of upperbound parameters values",
    )
    parser.add_argument(
        "-start",
        "--starttime",
        nargs="+",
        default=100,
        help="(int)\nStarting frame number for summary statistics",
    )
    parser.add_argument(
        "-end",
        "--endtime",
        nargs="+",
        default=200,
        help="(int)\nFinal frame number for summary statistics",
    )
    # summary statistics to calculate
    parser.add_argument(
        "-ssopts",
        "--summstatsopts",
        nargs="+",
        default=["A", "B", "C", "D","F" ,"H"],
        help="(list)\nSummary statistics to calculate (see allium/summstats.py for more details)",
    )
    parser.add_argument(
        "-state",
        "--state",
        type=str,
        default="confluent",
        help="(str)\nOption for calculation of posterior (though both sets of summstats are saved)",
    )
    parser.add_argument(
        "-framerate",
        "--framerate",
        type=str,
        default=1,
        help="(float)\nValue for the framerate conversion from output frame to [hours]",
    )
    # posterior options
    parser.add_argument(
        "--pcalc",
        dest="pcalc",
        help="(bool)\nCalculating posterior based on simulation run",
    )
    parser.add_argument(
        "-popt",
        "--posterioropt",
        default=["flow", "mix"],
        help="(list)\nList of posterior architectures to use",
    )
    
    parser.add_argument(
        "-maxN",
        "--maxN",
        type=int,
        default=5000,
        help="(int)\nNumber of cells to stop simulation",
    )

    parser.set_defaults(log=False, sim=False, test=False, pcalc=True)

    args = parser.parse_args()

    out = main(args)
