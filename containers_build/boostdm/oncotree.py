import os
import json
import functools
import operator

import pandas as pd

from boostdm.globals import COHORTS_PATH, ONCOTREE_PATH


def namespace(tree):
    pool = []
    for k, v in tree.items():
        pool += [k] + v
    return set(pool)


class Oncotree:

    def __init__(self):

        self.stats_cohorts = pd.read_csv(COHORTS_PATH, sep="\t")
        with open(ONCOTREE_PATH, 'r') as f:
            self.tree = json.load(f)
        self.ttypes = namespace(self.tree) - {'STOP'}

    def get_cohorts(self, ttype):
        """
        Given a ttype type ttype, retrieve the name of all cohorts that are associated with that ttype type
        :param ttype: string, name of ttype type
        :return: list of names of the cohorts
        """
        if ttype not in self.ttypes:
            raise Exception(f'ttype type {ttype} is not in oncotree namespace')
        cohorts = [
            tuple(x) for x in self.stats_cohorts[self.stats_cohorts["CANCER_TYPE"] == ttype][
                ["COHORT", "PLATFORM"]
            ].values]
        if ttype not in self.tree:  # ttype is a leaf, then return cohorts gathered from self.stats_cohorts
            return cohorts
        # ttype is a parent, then do a recursive call to get_cohorts the children nodes
        for child in self.tree[ttype]:
            cohorts += self.get_cohorts(child)
        return cohorts

    def get_ttypes(self, ttype):
        """
        Given a ttype, retrieve the name of all ttype types that are leaves in the oncotree
        :param ttype: string, name of ttype type
        :return: list of names of the ttype types
        """
        if ttype not in self.ttypes:
            raise Exception(f'ttype type {ttype} is not in oncotree namespace')
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
        :param ttype: string, name of ttype type
        :return: name of the parent
        """
        if (ttype not in self.ttypes) or (ttype == 'STOP'):
            return None
        
        for parent, children in self.tree.items():
            if ttype in children:
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


def get_leaves(k):
    oncotree = Oncotree()
    if k not in oncotree.tree:
        return [k]
    else:
        return functools.reduce(operator.concat, list(map(get_leaves, oncotree.tree[k])), [])


if __name__ == '__main__':

    tree = Oncotree()
