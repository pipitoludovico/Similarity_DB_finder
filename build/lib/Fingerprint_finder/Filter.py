class Filter:
    def __init__(self, df, idCol):
        print("Now dropping duplicates from the dataframe.")
        self.idCol = idCol
        self.df = df
        self.dfclean = self.df.drop_duplicates(subset=['CanonicalSmiles'])

    def getFiltered(self):
        return self.dfclean.loc[:, (self.idCol, 'CanonicalSmiles')]
