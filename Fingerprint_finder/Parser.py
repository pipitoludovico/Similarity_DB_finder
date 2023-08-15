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
        Column, Operator, Criterium = None, None, None
        allowedSelections = ('molecular_weight', 'formal_charge', 'similarity')
        ap = argparse.ArgumentParser()
        ap.add_argument("-fp", '--fingerprint', required=True,
                        help="specify -fp and put your .pdb file.")
        ap.add_argument('-db', '--database', required=True,
                        help=' use -db to input the SMILES database you want to explore')
        ap.add_argument('-so', '--sort', required=False, default="False",
                        help=' use -so True if you want your results to be sorted by ascending similarity, False for descending fashion')
        ap.add_argument('-se', '--separator', required=False, type=str, default=" ",
                        help=' use -se to define the separator of your database [Default = " "]')
        ap.add_argument("-k", '--kill', required=False, action='store_true', default=False,
                        help="Stop the current process.")
        ap.add_argument("-f", '--filter', nargs='*', required=False,
                        help='Choose between "molecular_weight", "formal_charge" or  "similarity" and add a criteria: eg. -f "molecular_weight >= 250" or -f "formal_charge <= 0"  or -f "similarity >= 70" to filter your results')
        ap.add_argument("-o", '--output', type=int, required=False, help='choose how many data you want to display in the results [Default = 100]')
        ap.add_argument("-sl", '--slice', nargs='*', required=False, help='choose a slice of your database to be used for processing. E.g. -sl 100:-1 [Default = 0 :-1]')
        args = ap.parse_args()

        if args.kill is True:
            from os import system, kill
            system('val=$(<.mypid ) && kill -9 $val')
            kill(getpid(), signal.SIGKILL)

        if str(args.sort).upper().startswith('T'):
            sort = True
        else:
            sort = False

        if 't' in str(args.separator):
            args.separator = "\t"

        output = 100

        if args.output:
            output = int(args.output)

        if args.slice is None:
            sliceStart, sliceEnd = 0, -1
        else:
            slicer = str(args.slice[0]).split(":")
            sliceStart, sliceEnd = int(slicer[0]), int(slicer[1])

        if args.filter is None:
            filters = ['similarity', ">=", "75"]
            Column, Operator, Criterium = filters[0], filters[1], float(filters[2])/100
        else:
            filterArgs = args.filter[0]
            filters = (filterArgs.split(" "))
            if filters[0] not in allowedSelections:
                raise ValueError("Please choose between the following: ", allowedSelections)
            else:
                if str(filters[0]).startswith('simi'):
                    Column, Operator, Criterium = filters[0], filters[1], float(filters[2])/100
                else:
                    Column, Operator, Criterium = filters[0], filters[1], float(filters[2])
        if args.fingerprint is not None and args.database is not None and args.separator is not None:
            return args.fingerprint, args.database, sort, args.separator, Column, Operator, Criterium, output, sliceStart, sliceEnd
