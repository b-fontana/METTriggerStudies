# coding: utf-8

_all_ = [ 'WriteDAGManager', 'dag_outputs' ]

import os
import re
import atexit # https://stackoverflow.com/questions/865115/how-do-i-correctly-clean-up-a-python-object
import sys
sys.path.append("..")

from utils import utils
from condor.job_writer import JobWriter

from luigi_conf.luigi_cfg import cfg
lcfg = cfg() #luigi configuration

@utils.set_pure_input_namespace
def dag_outputs(args):
    return JobWriter.define_dag_output( localdir=args.localdir,
                                        tag=args.tag,
                                        name='workflow' )

class WriteDAGManager:
    def __init__(self, localdir, tag, jobs, mode='long'):
        if mode not in ('short', 'long'):
            raise ValueError('Mode {} is not supported.'.format(mode))
        self.mode = mode
        
        out = dag_outputs( {'localdir': localdir, 'tag': tag} )
        self.this_file = open(out, 'w')
        atexit.register(self.cleanup)
        
        self.write_configuration()

        job_keys = { 'HistosData', 'HistosMC',
                     'CountsData', 'CountsMC',
                     'HaddHistoData', 'HaddHistoMC',
                     'HaddCountsData', 'HaddCountsMC',
                     'EffSF', 'EffSFAgg', 'Discr' }
        if mode == 'long': #add steps not strictly required
            job_keys.update(('Union', 'Closure'))
        assert set(jobs.keys()) == job_keys

        self.jobs = jobs
        self.define_all_job_names(self.jobs)

    def build_job_id(self, job_path):
        jp = os.path.dirname(job_path)
        regex = re.compile('.+/{}/(.+)'.format(lcfg.analysis_folders['subm']))
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
        # histos to hadd for data
        self.write_parent_child_hierarchy( parents=[x for x in self.jobs['HistosData']],
                                           childs=[self.jobs['HaddHistoData'][0]] )

        # histos to hadd for MC
        self.write_parent_child_hierarchy( parents=[x for x in self.jobs['HistosMC']],
                                           childs=[self.jobs['HaddHistoMC'][0]] )
        self.new_line()

        # hadd aggregation for Data
        self.write_parent_child_hierarchy( parents=[self.jobs['HaddHistoData'][0]],
                                           childs=[self.jobs['HaddHistoData'][1]] )

        # hadd aggregation for MC
        self.write_parent_child_hierarchy( parents=[self.jobs['HaddHistoMC'][0]],
                                           childs=[self.jobs['HaddHistoMC'][1]] )
        self.new_line()

        # counts to add for data
        self.write_parent_child_hierarchy( parents=[x for x in self.jobs['CountsData']],
                                           childs=[self.jobs['HaddCountsData'][0]] )

        # counts to add for MC
        self.write_parent_child_hierarchy( parents=[x for x in self.jobs['CountsMC']],
                                           childs=[self.jobs['HaddCountsMC'][0]] )
        self.new_line()

        # counts add aggregation for Data
        self.write_parent_child_hierarchy( parents=[self.jobs['HaddCountsData'][0]],
                                           childs=[self.jobs['HaddCountsData'][1]] )

        # counts add aggregation for MC
        self.write_parent_child_hierarchy( parents=[self.jobs['HaddCountsMC'][0]],
                                           childs=[self.jobs['HaddCountsMC'][1]] )
        self.new_line()

        # efficiencies/scale factors draw and saving
        self.write_parent_child_hierarchy( parents=[self.jobs['HaddHistoData'][1],
                                                    self.jobs['HaddHistoMC'][1]],
                                           childs=self.jobs['EffSF'] )

        # efficiencies/scale factors aggregate output files into one per channel
        self.write_parent_child_hierarchy( parents=self.jobs['EffSF'],
                                           childs=self.jobs['EffSFAgg'] )
        self.new_line()
        
        # variable discriminator
        self.write_parent_child_hierarchy( parents=self.jobs['EffSF'],
                                           childs=[x for x in self.jobs['Discr']] )

        if self.mode == 'long':
            self.new_line()
            
            # union weights calculator
            self.write_parent_child_hierarchy( parents=[x for x in self.jobs['Discr']],
                                               childs=[x for x in self.jobs['Union']] )
            self.new_line()

            # hadd union efficiencies (only MC)
            self.write_parent_child_hierarchy( parents=[x for x in self.jobs['Union']],
                                               childs=[x for x in self.jobs['Closure']] )

# condor_submit_dag -no_submit diamond.dag
# condor_submit diamond.dag.condor.sub
# https://htcondor.readthedocs.io/en/latest/users-manual/dagman-workflows.html#optimization-of-submission-time

