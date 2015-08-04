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
import parser

SCALE_TAG_LIST = ['+-AdvDif',\
                  # '+-InitFMM_Cheb', \
                  '+-SL', \
                  '+-FMM', \
                  '+-Merge' \
                      ]
PROFILE_TAG_LIST = ['+-EvalTree', \
                    '+-Evaluation'\
                    ]

def post_process_profile_node(node):
    if not node:
        return
    for profile_tag in PROFILE_TAG_LIST:
        if profile_tag in node.title:
            imb_val = "{0:>10.4f}".format( float(node.values[2])/float(node.values[0]) )
            node.values.append(imb_val)
    return node

def post_process_profile_data(mydoc, file_pp, PRINT_HEADER = True):
    """
    post processing of data
    Arguments:
    - `output`:
    - `file_pp`:
    """
    eval_tree_counter = 0
    for node in mydoc.node_list:
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

def post_process_scaling_data(mydoc, file_pp, PRINT_HEADER = True):
    """
    post processing of data
    Arguments:
    - `output`:
    - `file_pp`:
    """
    ppnode_title = mydoc.np
    ppnode_values = dict()

    for scale_tag in SCALE_TAG_LIST:
        for node in mydoc.node_list:
            if scale_tag in node.title:
                print node.title
                ppnode_values[node.title] = node.values['t_avg']
                if '+-SL' in node.title:
                    ppnode_values['f/s_total'] = node.values['f/s_total']
                    # print node.values['f/s_total']
                break
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
    # PARSE PROFILE OUTPUT
    # PRINT HEADER
    mydoc = parser.pdoc(output)
    pp_func(mydoc, file_pp, PRINT_HEADER);

if __name__ == '__main__':
    if len(sys.argv) >= 2:
        raw_dir_name = sys.argv[1]
    else:
        sys.exit()
    raw_file_names = list_raw_files(raw_dir_name)
    print raw_file_names
    file_pp = open(os.path.join(raw_dir_name, "sscal"+'.pp'), 'w')

    PRINT_HEADER = True
    for file_name in raw_file_names:
        f = open(file_name, 'r')
        post_process(f, file_pp, post_process_scaling_data, PRINT_HEADER)
        PRINT_HEADER = False
