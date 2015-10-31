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
from collections import OrderedDict

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

pattern_np_string = "-np\s(\d+\.?\d*)"
pattern_np = re.compile(pattern_np_string)

class pnode(object):
    """
    """

    def __init__(self, title, values):
        """
        """
        self.title = title
        self.values = values

    def print_me(self, file_prf):
        """
        print one node of profile tree
        Arguments:
        - `profile_node`:
        - `file_prf`:
        """
        string_format = "{:<50}".format(self.title)
        for key, val in self.values.iteritems():
            string_format += "{:>10}".format(val)
        string_format += "\n"
        file_prf.write(string_format)

class pdoc(object):
    """
    """

    def __init__(self, output):
        """
        """
        self.node_list = []
        self.header = []
        self.target_points_list = []
        self.np = 0
        for line in output:
            if line.startswith('TRG_CNT_IN_TOT:'):
                trg_cnt_match = pattern_trg_value.findall(line)
                self.target_points_list.append([int(trg_cnt) for trg_cnt in trg_cnt_match])
            if line.startswith('# CMD:'):
                np_match = pattern_np.findall(line)
                print np_match
                self.np = np_match[0]
            # CATCH PROFILE HEADER
            header = self.__parse_profile_header(line)
            if len(header):
                self.header = header
            # CATCH PROFILE NODE
            node_data = self.__parse_profile_node(line)
            if node_data:
                node = pnode(node_data[0], OrderedDict(zip(self.header, node_data[1])))
                self.node_list.append(node)

    def print_me(self, file_out):
        self.__print_profile_header(self.header, file_out)
        for node in self.node_list:
            node.print_me(file_out)


    def __parse_profile_header(self, line):
        """
        """
        header_match = []
        if PROFILE_TAG_HEADER in line:
            header_match = pattern_prof_header.findall(line)
        return header_match

    def __print_profile_header(self, header_list, file_pr):
        if len(header_list) == 0:
            return
        string_format = "{:<50}".format('# FUNCTION')
        for key, val in self.node_list[0].values.iteritems():
            string_format += "{:>10}".format(key)
        # for header in header_list:
        #     string_format += "{:>10}".format(header)
        string_format += "\n"
        file_pr.write('# ============================================================================================================================================\n')
        file_pr.write(string_format)
        file_pr.write('# ============================================================================================================================================\n')

    def __parse_profile_node(self, line):
        title_match = pattern_prof_title.findall(line)
        if len(title_match):
            value_match = pattern_prof_value.findall(line)
            return (title_match[0], value_match)
        return None
