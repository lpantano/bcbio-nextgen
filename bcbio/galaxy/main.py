def _graph_cmd(subparsers):
    parser = subparsers.add_parser("graph",
                                   help="Generate system graphs "
                                        "(CPU/memory/network/disk I/O "
                                        "consumption) from bcbio runs",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("log",
                        help="Local path to bcbio log file written by the run.")
    parser.add_argument("-o", "--outdir", default="monitoring/graphs",
                        help="Directory to write graphs to.")
    parser.add_argument("-r", "--rawdir", default="monitoring/collectl",
                        help="Directory to put raw collectl data files.")
    parser.add_argument("-v", "--verbose", action="store_true", default=False,
                        help="Emit verbose output")
    parser.set_defaults(func=graph.bootstrap)


