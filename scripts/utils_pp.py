#!/bin/env python

#*************************************************************************
#Copyright (C) 2015 by Arash Bakhtiari
#You may not use this file except in compliance with the License.
#You obtain a copy of the License in the LICENSE file.

#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.
#*************************************************************************
import re
import parser
import glob
import sys
import os
import time
import parser
import copy
from collections import OrderedDict
import json

def pp_profile_data(mydoc, file_pp, PRINT_HEADER = True):
    """
    post processing of data
    Arguments:
    - `output`:
    - `file_pp`:
    """

    PROFILE_TAG_LIST = [\
                    # '+-SL', \
                    # '+-EvalTree', \
                    '+-InEvaluation'\
                    ]

    ##############################
    # OUTOUT THE COMMAND
    ##############################
    file_pp.write(mydoc.cmd)
    print mydoc.cmd
    eval_tree_counter = 0
    ##############################
    # ITERATE OVER ALL PROFILE NODES
    ##############################
    for node in mydoc.node_list:
        # ITERATE OVER PROFILE TAGS
        for profile_tag in PROFILE_TAG_LIST:
            if profile_tag in node.title:
                time_imb = "{0:>10.2f}".format( float(node.values['t_max'])/float(node.values['t_min']) )
                node.values['t_imb'] = time_imb
                if '+-EvalTree' in node.title:
                    ##############################
                    # COMPUTE THE EVAL TREE MIN TIME DIFFERENTELY
                    ##############################
                    node_indx = mydoc.node_list.index(node)
                    t_min = 0.0;
                    for i in range(1, 6):
                        # print mydoc.node_list[node_indx+i].title
                        t_min += float(mydoc.node_list[node_indx+i].values['t_min'])
                    node.values['t_min'] = str(t_min)
                    ##############################
                    # CALCULATE THE TIME IMBALANCE
                    ##############################
                    time_imb = "{0:>10.2f}".format( float(node.values['t_max'])/float(node.values['t_min']))
                    node.values['t_imb'] = time_imb
                    ##############################
                    # CALCULATE THE PARTICLES IMBALANCE
                    ##############################
                    if len(mydoc.target_points_list) != 0:
                        max_points = max(mydoc.target_points_list[eval_tree_counter])
                        min_points = min(mydoc.target_points_list[eval_tree_counter])
                        target_points_imb = "{0:>10.2f}".format( float(max_points)/ min_points) 
                        node.values['p_imb'] = target_points_imb
                        eval_tree_counter += 1
                ##############################
                # PRINT HEADER/NODE
                ##############################
                if PRINT_HEADER:
                    header_string_format = "{:<50}".format('# FUNCTION')
                    for key, val in node.values.iteritems():
                        header_string_format += "{:>10}".format(key)
                    header_string_format += "\n"
                    file_pp.write(header_string_format)
                    PRINT_HEADER = False
                node.print_me(file_pp)

def pp_perf_cubic(mydoc, file_pp, PRINT_HEADER = True):
    """
    post processing of data
    Arguments:
    - `output`:
    - `file_pp`:
    """
    ##############################
    # print command
    ##############################
    print '--> post processing command:\n'+mydoc.cmd
    ##############################
    # add interpolation time average
    ##############################
    t_min_sum = 0
    t_max_sum = 0
    t_avg_sum = 0
    counter = 0
    for node in mydoc.node_list:
        # ITERATE OVER PROFILE TAGS
        if '+-InEvaluation' in node.title:
            counter = counter + 1
            t_min_sum += float(node.values['t_min'])
            t_avg_sum += float(node.values['t_avg'])
            t_max_sum += float(node.values['t_max'])
    ppnode_title  = '+-InEvaluation'
    ppnode_values = OrderedDict()
    if counter:
        ppnode_values['t_min'] = "{0:>10.2f}".format(t_min_sum/counter)
        # ppnode_values['t_avg'] = "{0:>10.2f}".format(t_avg_sum/counter)
        # ppnode_values['t_max'] = "{0:>10.2f}".format(t_max_sum/counter)
    ##############################
    # add the leaves count average
    ##############################
    # avoiding the float casting here on purpose
    if len(mydoc.leaves_count_list):
        ppnode_values['lvs_cnt'] = sum(mydoc.leaves_count_list)/len(mydoc.leaves_count_list)
    ##############################
    # create pp node to print
    ##############################
    ppnode = parser.pnode(ppnode_title, ppnode_values)
    ppnode.print_me(file_pp)

def pp_scal_s(mydoc, file_pp, PRINT_HEADER = True):
    """
    post processing of data
    Arguments:
    - `output`:
    - `file_pp`:
    """
    SCALE_TAG_LIST = [#'+-Solve',
                      '+-SLM',
                      '+-FMM',
                      'Merge',
                      '+-RefineTree',
                      '+-Balance21'
                      ]
    ppnode_title  = mydoc.np
    solve_tavg_sum = 0.0
    solve_favg_sum = 0.0
    tn_counter = 0

    ppnode_values = OrderedDict()
    ppnode_values['NP'] = mydoc.np
    for scale_tag in SCALE_TAG_LIST:
        ppnode_values['T'+scale_tag.replace('+-','')] = ''
        ppnode_values['F'+scale_tag.replace('+-','')] = ''
        for node in mydoc.node_list:
            if scale_tag in node.title:
                print node.title
                ppnode_values['T'+scale_tag.replace('+-','')] = node.values['t_avg']
                ppnode_values['F'+scale_tag.replace('+-','')] = node.values['f/s_total']
                solve_tavg_sum = solve_tavg_sum + float(node.values['t_avg'])
                solve_favg_sum = solve_favg_sum + float(node.values['f_avg'])
                tn_counter = tn_counter + 1
                break
    header_format = ''
    values_format = ''

    # ppnode_values['TSUM'] = "{:.4f}".format(solve_tavg_sum if tn_counter else 0)
    # ppnode_values['FSUM'] = "{:.4f}".format(solve_favg_sum if tn_counter else 0)
    # ppnode_values['FLOPS'] = "{:.4f}".format(solve_favg_sum/solve_tavg_sum  if tn_counter else 0)

    for key, val in ppnode_values.iteritems():
        header_format += "{:>12}".format(key)
        values_format += "{:>12}".format(val)
    header_format += "\n"
    values_format += "\n"
    file_pp.write(header_format)
    file_pp.write(values_format)


def pp_scal_s_mod(mydoc, file_pp, PRINT_HEADER = True):
    """
    post processing of data
    Arguments:
    - `output`:
    - `file_pp`:
    """
    # SCALE_TAG_LIST = [#'+-Solve',
    #                   '+-SLM',
    #                   '+-FMM',
    #                   'Merge',
    #                   '+-RefineTree',
    #                   '+-Balance21'
                      # ]

    # ppnode_title  = mydoc.np
    # solve_tavg_sum = 0.0
    # solve_favg_sum = 0.0
    # tn_counter = 0
    teval_sum = 0.0
    feval_sum = 0.0
    tsort_sum = 0.0
    tref_sum = 0.0
    tsolve_sum = 0.0
    fsolve_sum = 0.0
    tfmm = 0.0
    ffmm = 0.0

    ppnode_values = OrderedDict()
    ppnode_values['NP'] = mydoc.np
    for node in mydoc.node_list:
        scale_tag = '+-RefineTree'
        if  scale_tag in node.title:
            # ppnode_values['T'+scale_tag.replace('+-','')] = node.values['t_avg']
            tref_sum =  float(node.values['t_max'])

        scale_tag = 'Merge'
        if  scale_tag in node.title:
            # ppnode_values['T'+scale_tag.replace('+-','')] = node.values['t_max']
            tref_sum = tref_sum + float(node.values['t_max'])
            # tsolve_sum = tsolve_sum + float(node.values['t_max'])
            # fsolve_sum = fsolve_sum + float(node.values['f_avg'])

        scale_tag = 'Solve'
        if  scale_tag in node.title:
            # ppnode_values['T'+scale_tag.replace('+-','')] = node.values['t_max']
            tsolve_sum = tsolve_sum + float(node.values['t_max'])
            fsolve_sum = fsolve_sum + float(node.values['f/s_total'])

        # FMM
        scale_tag = '+-FMM'
        if  scale_tag in node.title:
            tfmm = float(node.values['t_max'])
            ffmm = float(node.values['f/s_total'])

        scale_tag = '+-InEvaluation'
        if  scale_tag in node.title:
            teval_sum = teval_sum + float(node.values['t_max'])
            feval_sum = feval_sum + float(node.values['f_avg'])

        scale_tag = '+-OutEvaluation'
        if  scale_tag in node.title:
            teval_sum = teval_sum + float(node.values['t_max'])
            feval_sum = feval_sum + float(node.values['f_avg'])

        scale_tag = '+-LclHQSort'
        if  scale_tag in node.title:
            tsort_sum = tsort_sum + float(node.values['t_max'])

        scale_tag = '+-OutScatterIndex'
        if  scale_tag in node.title:
            tsort_sum = tsort_sum + float(node.values['t_max'])

        scale_tag = '+-OutScatterForward'
        if  scale_tag in node.title:
            tsort_sum = tsort_sum + float(node.values['t_max'])

        scale_tag = '+-OutScatterReverse'
        if  scale_tag in node.title:
            tsort_sum = tsort_sum + float(node.values['t_max'])

    nn =    float( ppnode_values['NP'])
    ppnode_values['TSOLVE'] = "{:.4f}".format(tsolve_sum)
    ppnode_values['FSOLVE'] = "{:.4f}".format(fsolve_sum)
    ppnode_values['FSOLVEN'] = "{:.4f}".format(fsolve_sum/nn if nn else 0)
    ppnode_values['TREF']   = "{:.4f}".format(tref_sum)
    ppnode_values['TSORT']  = "{:.4f}".format(tsort_sum)
    ppnode_values['TEVAl']  = "{:.4f}".format(teval_sum)
    feval = feval_sum/teval_sum if teval_sum else 0
    ppnode_values['FEVAl']  = "{:.4f}".format(feval_sum/teval_sum if teval_sum else 0)
    ppnode_values['FEVAlN']  = "{:.4f}".format(feval/nn if nn else 0)
    ppnode_values['TFMM']   = "{:.4f}".format(tfmm)
    ppnode_values['FFMM']   = "{:.4f}".format(ffmm)
    ppnode_values['FFMMN']  = "{:.4f}".format(ffmm/nn if nn else 0)


    header_format = ''
    values_format = ''

    for key, val in ppnode_values.iteritems():
        header_format += "{:>12}".format(key)
        values_format += "{:>12}".format(val)
    header_format += "\n"
    values_format += "\n"
    file_pp.write(header_format)
    file_pp.write(values_format)

def pp_scal_s_mod_2(mydoc, file_pp, PRINT_HEADER = True):
    """
    post processing of data
    Arguments:
    - `output`:
    - `file_pp`:
    """
    # SCALE_TAG_LIST = [#'+-Solve',
    #                   '+-SLM',
    #                   '+-FMM',
    #                   'Merge',
    #                   '+-RefineTree',
    #                   '+-Balance21'
                      # ]

    # ppnode_title  = mydoc.np
    # solve_tavg_sum = 0.0
    # solve_favg_sum = 0.0
    # tn_counter = 0
    teval_sum = 0.0
    tmerge_sum = 0.0
    feval_sum = 0.0
    tsort_sum = 0.0
    tcomm_sum = 0.0
    tref_sum = 0.0
    tsolve_sum = 0.0
    fsolve_sum = 0.0
    tfmm = 0.0
    ffmm = 0.0
    tslm = 0.0
    fslm = 0.0

    ppnode_values = OrderedDict()
    ppnode_values['NP'] = mydoc.np
    for node in mydoc.node_list:
        scale_tag = '+-RefineTree'
        if  scale_tag in node.title:
            # ppnode_values['T'+scale_tag.replace('+-','')] = node.values['t_avg']
            tref_sum =  float(node.values['t_max'])

        scale_tag = 'Merge'
        if  scale_tag in node.title:
            # ppnode_values['T'+scale_tag.replace('+-','')] = node.values['t_max']
            tmerge_sum = tref_sum + float(node.values['t_max'])

        scale_tag = 'Solve'
        if  scale_tag in node.title:
            # ppnode_values['T'+scale_tag.replace('+-','')] = node.values['t_max']
            tsolve_sum = tsolve_sum + float(node.values['t_max'])
            fsolve_sum = fsolve_sum + float(node.values['f/s_total'])

        # FMM
        scale_tag = '+-FMM'
        if  scale_tag in node.title:
            tfmm = float(node.values['t_max'])
            ffmm = float(node.values['f/s_total'])

        # SLM
        scale_tag = '+-SLM'
        if  scale_tag in node.title:
            tslm = float(node.values['t_max'])
            fslm = float(node.values['f/s_total'])

        scale_tag = '+-InEvaluation'
        if  scale_tag in node.title:
            teval_sum = teval_sum + float(node.values['t_max'])
            feval_sum = feval_sum + float(node.values['f_avg'])
        scale_tag = '+-OutEvaluation'
        if  scale_tag in node.title:
            teval_sum = teval_sum + float(node.values['t_max'])
            feval_sum = feval_sum + float(node.values['f_avg'])

        scale_tag = '+-LclHQSort'
        if  scale_tag in node.title:
            tsort_sum = tsort_sum + float(node.values['t_max'])

        scale_tag = '+-OutScatterIndex'
        if  scale_tag in node.title:
            tcomm_sum = tcomm_sum + float(node.values['t_max'])
        scale_tag = '+-OutScatterForward'
        if  scale_tag in node.title:
            tcomm_sum = tcomm_sum + float(node.values['t_max'])
        scale_tag = '+-OutScatterReverse'
        if  scale_tag in node.title:
            tcomm_sum = tcomm_sum + float(node.values['t_max'])

    nn =    float( ppnode_values['NP'])
    ppnode_values['TSOLVE'] = "{:.4f}".format(tsolve_sum)
    ppnode_values['FSOLVE'] = "{:.4f}".format(fsolve_sum)
    ppnode_values['FSOLVEN'] = "{:.4f}".format(fsolve_sum/nn if nn else 0)

    ppnode_values['TREF']   = "{:.4f}".format(tref_sum)
    ppnode_values['TSORT']  = "{:.4f}".format(tsort_sum)

    ppnode_values['TEVAl']  = "{:.4f}".format(teval_sum)
    feval = feval_sum/teval_sum if teval_sum else 0
    ppnode_values['FEVAl']  = "{:.4f}".format(feval)
    # ppnode_values['FEVAlN']  = "{:.4f}".format(feval/nn if nn else 0)

    ppnode_values['TMRG']   = "{:.4f}".format(tmerge_sum)
    ppnode_values['TCOMM']  = "{:.4f}".format(tcomm_sum)

    ppnode_values['TFMM']   = "{:.4f}".format(tfmm)
    ppnode_values['FFMM']   = "{:.4f}".format(ffmm)
    ppnode_values['FFMMN']  = "{:.4f}".format(ffmm/nn if nn else 0)

    ppnode_values['TSLM']   = "{:.4f}".format(tslm)
    ppnode_values['FSLM']   = "{:.4f}".format(fslm)
    ppnode_values['FSLMN']  = "{:.4f}".format(fslm/nn if nn else 0)

    header_format = ''
    values_format = ''

    for key, val in ppnode_values.iteritems():
        header_format += "{:>12}".format(key)
        values_format += "{:>12}".format(val)
    header_format += "\n"
    values_format += "\n"
    file_pp.write(header_format)
    file_pp.write(values_format)


def get_time(mydoc):
    ppnode_values = OrderedDict()
    solve_tavg_sum = 0.0
    tn_counter = 0
    for node in mydoc.node_list:
        if "Solve_TN" in node.title:
            # print node.title
            solve_tavg_sum = solve_tavg_sum + float(node.values['t_avg'])
            tn_counter = tn_counter + 1
    ppnode_values['T_SUM'] = solve_tavg_sum if tn_counter else 0
    ppnode_values['T_AVG'] = solve_tavg_sum/tn_counter if tn_counter else 0
    return ppnode_values

def pp_scal_ws(mydoc, file_pp, PRINT_HEADER = True):
    """
    post processing of scal-ws.py outputs
    Arguments:
    - `output`:
    - `file_pp`:
    """
    ##############################
    # OUTOUT THE COMMAND
    ##############################
    print '--> post processing command:\n'+mydoc.cmd

    ##############################
    # CREATE THE PP LINE
    ##############################
    ppnode_title  = mydoc.np
    ppnode_values = get_time(mydoc)

    ##############################
    # FORMAT/PRINT THE PP LINE
    ##############################
    ppnode = parser.pnode(ppnode_title, ppnode_values)
    titlew = 15
    valuew = 15
    title_format = "{:>"+str(titlew)+"}"
    title_value_format = "{:>"+str(valuew)+"}"
    value_format = "{:>"+str(valuew)+".2f}"
    if PRINT_HEADER:
        header_string_format = title_format.format('NP')
        for key, val in ppnode_values.iteritems():
            header_string_format += title_value_format.format(key)
        header_string_format += "\n"
        file_pp.write(header_string_format)
    ppnode.print_me(file_pp, title_format, value_format)
    file_pp.flush()

def pp_tree_eval_data(mydoc, file_pp, PRINT_HEADER = True):
    """
    post processing of data
    Arguments:
    - `output`:
    - `file_pp`:
    """
    SCALE_TAG_LIST = ['+-LclSort',\
                      '+-GlobalSort'
                      ]

    ppnode_title = mydoc.np
    ppnode_values = OrderedDict()

    for scale_tag in SCALE_TAG_LIST:
        for node in mydoc.node_list:
            if scale_tag in node.title:
                print node.title
                ppnode_values[node.title] = node.values['t_avg']
                # if '+-SL' in node.title:
                #     ppnode_values['f/s_total'] = node.values['f/s_total']
                #     print node.values['f/s_total']
                # break
    ppnode = parser.pnode(ppnode_title, ppnode_values)
    if PRINT_HEADER:
        header_string_format = "{:<50}".format('NP')
        for key, val in ppnode_values.iteritems():
            header_string_format += "{:>10}".format(key)
        header_string_format += "\n"
        file_pp.write(header_string_format)
    ppnode.print_me(file_pp)


def list_raw_files(raw_dir_name, extension='.out'):
    """
    """
    import os
    outfiles = [os.path.join(root, filename) \
              for root, dirnames, filenames in os.walk(raw_dir_name) \
                for filename in filenames if filename.endswith(extension)]
    return sorted(outfiles)
    # return sorted(glob.glob(raw_dir_name + "*.out"))

def post_process(raw_files_dir, pp_func):
    ############################################################
    # COLLECT A LIST OF ALL AVAIABLE RAW OUTPUT FILES
    ############################################################
    raw_files_list = list_raw_files(raw_files_dir)
    ############################################################
    # CREATE/OPEN THE PP OUTPUT FILE
    ############################################################
    TIMESTR       = time.strftime("%Y%m%d-%H%M%S")
    pp_output_path = os.path.join(raw_files_dir, 'pp_output_'+pp_func.__name__+'_'+TIMESTR+'.pp')
    fpp = open(pp_output_path, 'w')
    ############################################################
    # POST PROCESS FILES
    ############################################################
    # print '--> printing commands'
    # [fpp.write(parser.pdoc(open(raw_file, 'r')).cmd)  for raw_file in raw_files_list]
    PRINT_HEADER = True
    for raw_file in raw_files_list:
        pp_func(parser.pdoc(open(raw_file, 'r')), fpp, PRINT_HEADER)
        # PRINT_HEADER = False
    # [pp_func(parser.pdoc(open(raw_file, 'r')), fpp, PRINT_HEADER) for raw_file in raw_files_list]

if __name__ == '__main__':
    if len(sys.argv) >= 2:
        raw_files_dir = sys.argv[1]
    else:
        sys.exit()

    # pp_func = pp_scal_s
    pp_func = pp_scal_s_mod_2
    # pp_func = pp_tree_eval_data
    # pp_func = pp_profile_data
    # pp_func = pp_perf_cubic
    # pp_func = pp_scal_ws
    post_process(raw_files_dir, pp_func)
