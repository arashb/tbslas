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

PROFILE_TAG_HEADER = 't_min'
# PROFILE_TAG_LIST = ['+-RunSemilag', '+-RunFMM', '+-EvalTree', \
#                     '+-MortonId', '+-ScatterIndex', '+-ScatterForward',\
#                     '+-Evaluation', '+-ScatterReverse']
PROFILE_TAG_LIST = ['+-EvalTree', \
                    '+-Evaluation'\
                    ]

# POINT TARGET VALUES RE PATTERN
pattern_trg_value_string = "(\d+)"
pattern_trg_value = re.compile(pattern_trg_value_string)

# PROFILE HEADER RE PATTERN
pattern_prof_header_string = "(\w+[//]*\w+)"
pattern_prof_header = re.compile(pattern_prof_header_string)

# PROFILE VALUES RE PATTERNS
pattern_prof_value_string = "(\d+\.\d+)"
pattern_prof_value = re.compile(pattern_prof_value_string)

pattern_prof_title_string = "([\|\s]*\+\-\w+)"
pattern_prof_title = re.compile(pattern_prof_title_string)

class pdoc(object):
    """
    """

    def __init__(self, output):
        """
        """
        self.node_list = []
        self.header = []
        self.target_points_list = []
        for line in output:
            if line.startswith('TRG_CNT'):
                trg_cnt_match = pattern_trg_value.findall(line)
                self.target_points_list.append([int(trg_cnt) for trg_cnt in trg_cnt_match])

            header = parse_profile_header(line)
            if len(header):
                self.header = header
            # CATCH PROFILE NODE
            node = parse_profile_node(line)
            if node:
                self.node_list.append(node)

    def print_me(self, file_out):
        print_profile_header(self.header, file_out)
        for node in self.node_list:
            print_profile_node(node, file_out)

class pnode(object):
    """
    """

    def __init__(self, title, values):
        """
        """
        self.title = title
        self.values = values
        self.children = []

def parse_profile_header(line):
    """
    """
    header_match = []
    if PROFILE_TAG_HEADER in line:
            header_match = pattern_prof_header.findall(line)
    return header_match

def print_profile_header(header_list, file_pr):
    if len(header_list) == 0:
        return
    string_format = "{:<50}".format('# FUNCTION')
    for header in header_list:
        string_format += "{:>10}".format(header)
    string_format += "\n"
    file_pr.write('# ============================================================================================================================================\n')
    file_pr.write(string_format)
    file_pr.write('# ============================================================================================================================================\n')

def parse_profile_node(line):
    title_match = pattern_prof_title.findall(line)
    if len(title_match):
        value_match = pattern_prof_value.findall(line)
        return pnode(title_match[0], value_match)
    return None

def print_profile_node(node, file_prf):
    """
    print one node of profile tree
    Arguments:
    - `profile_node`:
    - `file_prf`:
    """
    if not node:
        return
    string_format = "{:<50}".format(node.title)
    for val in node.values:
        string_format += "{:>10}".format(val)
    string_format += "\n"
    file_prf.write(string_format)

def parse_profile_data(output, file_pr, PRINT_PRFL_HEADER):
    """
    catch and save profiling data
    Arguments:
    - `output`:
    - `file_pr`:
    """
    node_list = []
    for line in output:
        if PRINT_PRFL_HEADER:
            header_match = parse_profile_header(line)
            print_profile_header(header_match, file_pr)
        # CATCH PROFILE DATA
        node = parse_profile_node(line)
        node_list.append(node)
        print_profile_node(node, file_pr)
    return node_list
