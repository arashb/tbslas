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

pattern_np_string = "-np?\s(\d+\.?\d*)"
pattern_np = re.compile(pattern_np_string)

class pnode(object):
    """
    """

    def __init__(self, title, values):
        """
        """
        self.title = title
        self.values = values

    def print_me(self, file_prf, title_format="{:<50}", value_format="{:<10}"):
        """
        print one node of profile tree
        Arguments:
        - `profile_node`:
        - `file_prf`:
        """
        line_string = title_format.format(self.title)
        for key, val in self.values.iteritems():
            line_string += value_format.format(val)
        line_string += "\n"
        file_prf.write(line_string)

class pdoc(object):
    """
    """

    def __init__(self, output):
        """
        """
        self.node_list          = []
        self.header             = []
        self.target_points_list = []
        self.np                 = 0
        self.cmd                = ''
        self.reporter_header    = []
        self.reporter_results   = []
        self.reporter           = OrderedDict()
        self.leaves_count_list  = []
        self.max_mem_per_node   = None

        for line in output:
            # CAPTURE TOTAL TARGET POINTS COUNT
            if line.startswith('TRG_CNT_IN_TOT:'):
                trg_cnt_match = pattern_trg_value.findall(line)
                self.target_points_list.append([int(trg_cnt) for trg_cnt in trg_cnt_match])

            # CAPTURE THE RUNNING COMMAND + NUM PROCS
            if line.startswith('# CMD '):
                self.cmd = line
                np_match = pattern_np.findall(line)
                self.np = np_match[0]

            # CAPTURE TBSLAS HEADER
            if line.startswith('#TBSLAS-HEADER:'):
                self.reporter_header = line.split()[1:]    # remove the line tag (first element)

            # CAPTURE TBSLAS RESULT
            if line.startswith('#TBSLAS-RESULT:'):
                self.reporter_results = line.split()[1:]   # remove the line tag (first element)

            # CAPTURE LEAVES COUNT
            if line.startswith('# TOT_LEAVES_CNT'):
                self.leaves_count_list.append(int(line.split()[2]))
                # print line.split()[2:]

            if line.startswith("TACC: Max Memory Used Per Node"):
                self.max_mem_per_node = line.split(": ")[2]
            # CAPTURE PROFILE HEADER
            header = self.__parse_profile_header(line)
            if len(header):
                self.header = header

            # CAPTURE PROFILE NODE DATA
            node_data = self.__parse_profile_node(line)
            if node_data:
                node = pnode(node_data[0], OrderedDict(zip(self.header, node_data[1])))
                self.node_list.append(node)

        self.reporter = OrderedDict(zip(self.reporter_header, self.reporter_results))
        # print self.reporter
        # print self.leaves_count_list


    def print_me(self, file_out):
        self.__print_profile_header(self.header, file_out)
        for node in self.node_list:
            node.print_me(file_out)

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

    def __parse_profile_header(self, line):
        """
        """
        header_match = []
        if PROFILE_TAG_HEADER in line:
            header_match = pattern_prof_header.findall(line)
        return header_match

    def __parse_profile_node(self, line):
        title_match = pattern_prof_title.findall(line)
        if len(title_match):
            value_match = pattern_prof_value.findall(line)
            return (title_match[0], value_match)
        return None
