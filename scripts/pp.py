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


SCALE_TAG_LIST = ['+-AdvDif',\
                  # '+-InitFMM_Cheb', \
                  '+-SL', \
                  '+-FMM', \
                  '+-CMerge' \
                  '+-SMerge' \
                      ]

PROFILE_TAG_LIST = [\
                    # '+-SL', \
                    # '+-EvalTree', \
                    '+-InEvaluation'\
                    ]

def pp_profile_node(node):
    if not node:
        return
    for profile_tag in PROFILE_TAG_LIST:
        if profile_tag in node.title:
            imb_val = "{0:>10.4f}".format( float(node.values[2])/float(node.values[0]) )
            node.values.append(imb_val)
    return node

def pp_profile_data(mydoc, file_pp, PRINT_HEADER = True):
    """
    post processing of data
    Arguments:
    - `output`:
    - `file_pp`:
    """
    # OUTOUT THE COMMAND
    file_pp.write(mydoc.cmd)
    print mydoc.cmd

    eval_tree_counter = 0
    # ITERATE OVER ALL PROFILE NODES
    for node in mydoc.node_list:
        # ITERATE OVER PROFILE TAGS
        for profile_tag in PROFILE_TAG_LIST:
            if profile_tag in node.title:
                time_imb = "{0:>10.4f}".format( float(node.values['t_max'])/float(node.values['t_min']) )
                node.values['t_imb'] = time_imb
                if '+-EvalTree' in node.title:
                    # COMPUTE THE EVAL TREE MIN DIFFERENTELY
                    node_indx = mydoc.node_list.index(node)
                    t_min = 0.0;
                    for i in range(1, 6):
                        # print mydoc.node_list[node_indx+i].title
                        t_min += float(mydoc.node_list[node_indx+i].values['t_min'])
                    node.values['t_min'] = str(t_min)
                    time_imb = "{0:>10.4f}".format( float(node.values['t_max'])/float(node.values['t_min']))
                    node.values['t_imb'] = time_imb
                    if len(mydoc.target_points_list) != 0:
                        max_points = max(mydoc.target_points_list[eval_tree_counter])
                        min_points = min(mydoc.target_points_list[eval_tree_counter])
                        target_points_imb = "{0:>10.4f}".format( float(max_points)/ min_points) 
                        node.values['p_imb'] = target_points_imb
                        eval_tree_counter += 1
                if PRINT_HEADER:
                    header_string_format = "{:<50}".format('# FUNCTION')
                    for key, val in node.values.iteritems():
                        header_string_format += "{:>10}".format(key)
                    header_string_format += "\n"
                    file_pp.write(header_string_format)
                    PRINT_HEADER = False
                node.print_me(file_pp)

def pp_cubic_perf_data(mydoc, file_pp, PRINT_HEADER = True):
    """
    post processing of data
    Arguments:
    - `output`:
    - `file_pp`:
    """
    ##############################
    # print command
    ##############################
    file_pp.write(mydoc.cmd)
    print mydoc.cmd

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
        ppnode_values['t_min'] = "{0:>10.4f}".format(t_min_sum/counter)
        # ppnode_values['t_avg'] = "{0:>10.4f}".format(t_avg_sum/counter)
        # ppnode_values['t_max'] = "{0:>10.4f}".format(t_max_sum/counter)

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

def pp_scaling_data(mydoc, file_pp, PRINT_HEADER = True):
    """
    post processing of data
    Arguments:
    - `output`:
    - `file_pp`:
    """
    ppnode_title  = mydoc.np
    ppnode_values = OrderedDict()

    for scale_tag in SCALE_TAG_LIST:
        for node in mydoc.node_list:
            if scale_tag in node.title:
                print node.title
                ppnode_values[node.title] = node.values['t_avg']
                if '+-SL' in node.title:
                    ppnode_values['f/s_total'] = node.values['f/s_total']
                    print node.values['f/s_total']
                break
    ppnode = parser.pnode(ppnode_title, ppnode_values)
    if PRINT_HEADER:
        header_string_format = "{:<50}".format('NP')
        for key, val in ppnode_values.iteritems():
            header_string_format += "{:>10}".format(key)
        header_string_format += "\n"
        file_pp.write(header_string_format)
    ppnode.print_me(file_pp)

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


def list_raw_files(raw_dir_name):
    """
    """
    import os
    outfiles = [os.path.join(root, filename) \
              for root, dirnames, filenames in os.walk(raw_dir_name) \
              for filename in filenames if filename.endswith('.out')]
    return sorted(outfiles)
    # return sorted(glob.glob(raw_dir_name + "*.out"))

def post_process(output, file_pp, pp_func, PRINT_HEADER):
    mydoc = parser.pdoc(output)
    pp_func(mydoc, file_pp, PRINT_HEADER);

if __name__ == '__main__':
    if len(sys.argv) >= 2:
        raw_files_dir = sys.argv[1]
    else:
        sys.exit()

    ############################################################
    # collect a list of all avaiable raw output files
    ############################################################
    raw_files_list = list_raw_files(raw_files_dir)
    print raw_files_list

    # create post processing output file
    TIMESTR       = time.strftime("%Y%m%d-%H%M%S")
    pp_output_path = os.path.join(raw_files_dir, 'pp_output_'+TIMESTR+'.pp')
    pp_output_path = os.path.join(raw_files_dir, 'pp_output'+'.pp')

    pp_output_file = open(pp_output_path, 'w')

    PRINT_HEADER = True
    for raw_file in raw_files_list:
        f = open(raw_file, 'r')
        # post_process(f, pp_output_file, pp_scaling_data, PRINT_HEADER)
        # post_process(f, pp_output_file, pp_tree_eval_data, PRINT_HEADER)
        # post_process(f, pp_output_file, pp_profile_data, PRINT_HEADER)
        post_process(f, pp_output_file, pp_cubic_perf_data, PRINT_HEADER)
        PRINT_HEADER = False

    pp_output_file.close()
