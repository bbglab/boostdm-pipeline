import os
import json

from boostdm.globals import COHORTS_PATH, ONCOTREE_PATH

import pandas as pd


def namespace(tree):
    pool = []
    for k, v in tree.items():
        pool += [k] + v
    return set(pool)


class Oncotree:

    def __init__(self):

        self.stats_cohorts = pd.read_csv(COHORTS_PATH, sep='\t')
        self.tree_table = pd.read_csv(ONCOTREE_PATH, sep='\t')
        self.tree = self.to_dict()
        self.ttypes = namespace(self.tree)


    def to_dict(self):
        d = {}
        children_parent_pairs = list(zip(self.tree_table['ID'], self.tree_table['PARENT']))
        for child, parent in children_parent_pairs:
            if child != 'CANCER':
                d[parent] = d.get(parent, []) + [child]
        return d


    def get_cohorts(self, ttype):
        """
        Given a tumor type ttype, retrieve the name of all cohorts that are associated with that tumor type
        :param ttype: string, name of tumor type
        :return: cohort list
        """
        if ttype not in self.ttypes:
            raise Exception(f'tumor type {ttype} is not in oncotree namespace')

        cohorts = self.stats_cohorts[self.stats_cohorts["CANCER_TYPE"] == ttype]['COHORT'].tolist()

        # cohorts = [
        # tuple(x) for x in self.stats_cohorts[self.stats_cohorts["CANCER_TYPE"] == ttype][
        #        #["COHORT", "SOURCE"]
        #        ['COHORT']
        #    ].values]

        if ttype not in self.tree:  # ttype is a leaf, then return cohorts gathered from self.stats_cohorts
            return cohorts

        # ttype is a parent, then do a recursive call to get_cohorts the children nodes
        for child in self.tree[ttype]:
            cohorts += self.get_cohorts(child)

        return cohorts

    
    def get_ttypes(self, ttype):
        """
        Given a tumor type ttype, retrieve the name of all tumor types that are leaves in the oncotree
        :param ttype: string, name of tumor type
        :return: list of names of the tumor types
        """
        if ttype not in self.ttypes:
            raise Exception(f'tumor type {ttype} is not in oncotree namespace')
        ttypes, res = [ttype], []
        while ttypes:
            tt = ttypes.pop(0)
            children = self.tree.get(tt, [])
            if len(children) > 0:
                ttypes += children
            else:
                res.append(tt)
        return res


    def fetch_parent_ttype(self, ttype):
        """
        For a given ttype retrieve the parent ttype
        :param ttype: string, name of tumor type
        :return: name of the parent
        """

        if (ttype not in self.ttypes) or (ttype == 'CANCER'):
            return None

        for parent, childs in self.tree.items():
            if ttype in childs:
                return parent

        return ttype


    def fetch_parent_cohort(self, cohort):
        """
        For a given cohort retrieve the parent ttype
        :param cohort: string, name of COHORT
        :return: name of the parent
        """
        parent = self.stats_cohorts[self.stats_cohorts["COHORT"] == cohort]["CANCER_TYPE"].unique()
        if len(parent) > 0:
            return parent[0]
        else:
            return None


    def get_parents(self):
        return self.tree.keys()


    def is_cohort(self, cohort):
        return cohort in self.stats_cohorts["COHORT"].unique()


if __name__ == '__main__':

    tree = Oncotree()
    print(tree.fetch_parent_ttype("LGG"))
    print(tree.get_cohorts("CANCER"))
    print(tree.fetch_parent_cohort("ICGC_WXS_BOCA_UK"))
