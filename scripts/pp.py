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

# PROFILE_TAG_LIST = ['+-RunSemilag', '+-RunFMM', '+-EvalTree', \
#                     '+-MortonId', '+-ScatterIndex', '+-ScatterForward',\
#                     '+-Evaluation', '+-ScatterReverse']
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

def post_process_profile_data(mydoc, file_pp, PRINT_PRFL_HEADER = True):
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
                time_imb = "{0:>10.4f}".format( float(node.values[2])/float(node.values[0]) )
                node.values.append(time_imb)
                if '+-EvalTree' in node.title:
                    max_points = max(mydoc.target_points_list[eval_tree_counter])
                    min_points = min(mydoc.target_points_list[eval_tree_counter])
                    target_points_imb = "{0:>10.4f}".format( float(max_points)/ min_points) 
                    node.values.append(target_points_imb)
                    eval_tree_counter += 1
                parser.print_profile_node(node, file_pp)
