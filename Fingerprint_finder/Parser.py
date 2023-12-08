import argparse
from os import getpid
import signal


class ArgParser:
    def __init__(self):
        self.GetPID()
        self.ParseArgs()

    @staticmethod
    def GetPID():
        _ = getpid()
        pidFile = open(".mypid", "w")
        pidFile.write(str(_))
        pidFile.close()

    @staticmethod
    def ParseArgs():
        allowedSelections = ('molecular_weight', 'formal_charge', 'similarity')
        ap = argparse.ArgumentParser()
        ap.add_argument("-fp", '--fingerprint', required=True,
                        help="specify -fp and put your .pdb file.")
        ap.add_argument('-db', '--database', required=True,
                        help=' use -db to input the SMILES database you want to explore')
        ap.add_argument('-so', '--sort', required=False, default="False",
                        help=' use -so True if you want your results to be sorted by ascending similarity, False for descending fashion')
        ap.add_argument("-k", '--kill', required=False, action='store_true', default=False,
                        help="Stop the current process.")
        ap.add_argument('-f', '--filter', action='append', required=False, help='Set filters for your dataframe between between: molecular_weight, formal_charge, similarity. E.g.: -f "molecular_weight <= 455" -f "formal_charge <= 0"...')
        ap.add_argument("-o", '--output', type=int, required=False,
                        help='choose how many data you want to display in the results [Default = 100]')
        ap.add_argument("-sl", '--slice', nargs='*', required=False,
                        help='choose a slice of your database to be used for processing. E.g. -sl 100:-1 [Default = 0 :-1]')
        ap.add_argument('-ex', '--exclude', action='append', required=False, help='Exclude all rows containing the described smiles. E.g.: -ex "O-" -ex "S(=O)"')
        ap.add_argument('-in', '--include', action='append', required=False, help='Include only the rows containing the described smiles. E.g.: -in "O-" -in "S(=O)"')
        args = ap.parse_args()

        if args.kill is True:
            from os import system, kill
            ArgParser().GetPID()
            system('val=$(<.mypid ) && kill -9 $val')
            kill(getpid(), signal.SIGKILL)

        if str(args.sort).upper().startswith('T'):
            sort = True
        else:
            sort = False

        output = 100

        if args.output:
            output = int(args.output)

        if args.slice is None:
            sliceStart, sliceEnd = 0, 100
        else:
            slicer = str(args.slice[0]).split(":")
            sliceStart, sliceEnd = int(slicer[0]), int(slicer[1])

        if args.filter is None:
            args.filter = ["similarity >= 75"]

        return args.fingerprint, args.database, sort, args.filter, output, sliceStart, sliceEnd, args.exclude, args.include
