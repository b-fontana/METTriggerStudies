# coding: utf-8

_all_ = [ 'WriteDAGManager', 'dag_outputs' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

import re
import atexit # https://stackoverflow.com/questions/865115/how-do-i-correctly-clean-up-a-python-object

import inclusion
from inclusion.utils import utils
from inclusion.condor.job_writer import JobWriter
from inclusion.config import main

@utils.set_pure_input_namespace
def dag_outputs(args):
    return JobWriter.define_dag_output( localdir=args.localdir,
                                        tag=args.tag,
                                        name='workflow' )

class WriteDAGManager:
    def __init__(self, localdir, tag, jobs, branch='all'):
        if branch not in ('all', 'counts', 'extra'):
            raise ValueError('Branch {} is not supported.'.format(branch))
        self.branch = branch
        
        out = dag_outputs( {'localdir': localdir, 'tag': tag} )
        self.this_file = open(out, 'w')
        atexit.register(self.cleanup)
        
        self.write_configuration()

        job_keys = {'HistosData', 'HistosMC',
                    'CountsData', 'CountsMC',
                    'HaddHistoData', 'HaddHistoMC',
                    'HaddCountsData', 'HaddCountsMC',
                    'EffSF', 'EffSFAgg', 'Discr',
                    'Union', 'Closure'}
        for k in jobs:
            assert k in job_keys

        self.jobs = jobs
        self.define_all_job_names(self.jobs)

    def build_job_id(self, job_path):
        jp = os.path.dirname(job_path)
        regex = re.compile('.+/{}/(.+)'.format(main.folders['subm']))
        matches = regex.findall(jp)
        assert len(matches) == 1
        matches = matches[0].replace('/', '_')
        return matches
         
    def cleanup(self):
        self.this_file.close()
            
    def write_string(self, string):
        self.this_file.write(string)
            
    def new_line(self):
        self.this_file.write('\n')

    def define_all_job_names(self, jobs):
        for _,values in jobs.items():
            self.define_job_names(values)
      
    def define_job_names(self, jobs):
        """First step to build a DAG"""
        for job in jobs:
            job_id = self.build_job_id(job)
            self.write_string('JOB {} {}\n'.format(job_id, job))
        self.new_line()

    def write_parent_child_hierarchy(self, parents, childs):
        if not isinstance(parents, (list,tuple)):
            m = ' Please pass lists to the '
            m += ' `write_parent_child_hierarchy method.'
            raise TypeError(m)
        
        self.write_string('PARENT ')
        for par in parents:
            job_id = self.build_job_id(par)
            self.write_string('{} '.format(job_id))
        self.write_string('CHILD ')
        for cld in childs:
            job_id = self.build_job_id(cld)
            self.write_string('{} '.format(job_id))
        self.new_line()

    def write_configuration(self):
        # https://research.cs.wisc.edu/htcondor/manual/v7.6/2_10DAGMan_Applications.html#SECTION003109000000000000000
        # dot -Tps dag.dot -o dag.ps
        self.this_file.write('DOT dag.dot\n')
        #self.this_file.write('DAGMAN_HOLD_CLAIM_TIME=30\n')
        self.new_line()

    def write_all(self):
        # counts to add for data
        p = [x for x in self.jobs['CountsData']]
        c = [self.jobs['HaddCountsData'][0]]
        self.write_parent_child_hierarchy(parents=p, childs=c)

        # counts to add for MC
        p = [x for x in self.jobs['CountsMC']]
        c = [self.jobs['HaddCountsMC'][0]]
        self.write_parent_child_hierarchy(parents=p, childs=c)
        self.new_line()

        # counts add aggregation for Data
        p = [self.jobs['HaddCountsData'][0]]
        c = [self.jobs['HaddCountsData'][1]]
        self.write_parent_child_hierarchy(parents=p, childs=c)

        # counts add aggregation for MC
        p = [self.jobs['HaddCountsMC'][0]]
        c = [self.jobs['HaddCountsMC'][1]]
        self.write_parent_child_hierarchy(parents=p, childs=c)
        self.new_line()

        if self.branch == 'all':
            # histos to hadd for dat
            p = [x for x in self.jobs['HistosData']]
            c = [self.jobs['HaddHistoData'][0]]
            self.write_parent_child_hierarchy(parents=p, childs=c)

            # histos to hadd for MC
            p = [x for x in self.jobs['HistosMC']]
            c = [self.jobs['HaddHistoMC'][0]]
            self.write_parent_child_hierarchy(parents=p, childs=c)
            self.new_line()

            # hadd aggregation for Data
            p = [self.jobs['HaddHistoData'][0]]
            c = [self.jobs['HaddHistoData'][1]]
            self.write_parent_child_hierarchy(parents=p, childs=c)

            # hadd aggregation for MC
            p = [self.jobs['HaddHistoMC'][0]]
            c = [self.jobs['HaddHistoMC'][1]]
            self.write_parent_child_hierarchy(parents=p, childs=c)
            self.new_line()

            # efficiencies/scale factors draw and saving
            p = [self.jobs['HaddHistoData'][1], self.jobs['HaddHistoMC'][1]]
            c = self.jobs['EffSF']
            self.write_parent_child_hierarchy(parents=p, childs=c)

            # efficiencies/scale factors aggregate output files into one per channel
            p = self.jobs['EffSF']
            c = self.jobs['EffSFAgg']
            self.write_parent_child_hierarchy(parents=p, childs=c)
        
        if self.branch == 'extra':
            self.new_line()

            # variable discriminator
            p = self.jobs['EffSF']
            c = [x for x in self.jobs['Discr']]
            self.write_parent_child_hierarchy(parents=p, childs=c)
            self.new_line()
                
            # union weights calculator
            p = [x for x in self.jobs['Discr']]
            c = [x for x in self.jobs['Union']]
            self.write_parent_child_hierarchy(parents=p, childs=c)
            self.new_line()

            # hadd union efficiencies (only MC)
            p = [x for x in self.jobs['Union']]
            c = [x for x in self.jobs['Closure']]
            self.write_parent_child_hierarchy(parents=p, childs=c)
