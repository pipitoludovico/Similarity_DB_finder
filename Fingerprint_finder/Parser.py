import argparse
from os import getpid
import signal


class ArgParser:
    def __init__(self):
        getpid()
        self.ParseArgs()

    @staticmethod
    def GetPID():
        _ = getpid()
        pidFile = open(".mypid", "w")
        pidFile.write(str(_))
        pidFile.close()

    @staticmethod
    def ParseArgs():
        ap = argparse.ArgumentParser()
        ap.add_argument("-fp", '--fingerprint', required=True,
                        help="specify -fp and put your .pdb file.")
        ap.add_argument('-db', '--database', required=True,
                        help=' use -db to input the SMILES database you want to explore')
        ap.add_argument("-k", '--kill', required=False, action='store_true', default=False,
                        help="Stop the current process.")

        args = ap.parse_args()

        if args.kill is True:
            from os import system, kill
            system('val=$(<.mypid ) && kill -9 $val')
            kill(getpid(), signal.SIGKILL)

        if args.fingerprint is not None and args.database is not None:
            return args.fingerprint, args.database