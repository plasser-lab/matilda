"""
Utilities for reading and writing options from/to an input file.
"""

from __future__ import print_function, division
import sys
from . import error_handler


def user_input_py2(inpstr):
    return raw_input(inpstr)


def user_input_py3(inpstr):
    return input(inpstr)


if sys.version_info[0] == 2:
    user_input = user_input_py2
else:
    user_input = user_input_py3


class options:
    """
    Base class for handling input options.
    """
    def __init__(self, ifile):
        self.opt_dict = {}
        self.descr_dict = {} # Collect information for documentation
        self.doc_list = []
        self.ifile = ifile

    def __getitem__(self, option):
        return self.get(option, strict=True)

    def get(self, option, strict=True):
        """
        Return the value of an option.
        """
        self.chk_option(option)

        if strict and self.opt_dict[option] == None:
            raise error_handler.MsgError('Option "%s" not defined in file %s!'%(option, self.ifile))
        else:
            return self.opt_dict[option]

    def __setitem__(self, key, val):
        self.opt_dict[key] = val

    def set_kd(self, key, val, descr=''):
        """
        Set key, value along with a description.
        """
        self.opt_dict[key]   = val
        self.descr_dict[key] = descr
        self.doc_list.append(key)

    def __contains__(self, option):
        """
        Check if an option has been set.

        Raise an error if the option does not even exist.
        """
        self.chk_option(option)

        return self.opt_dict[option] != None

    def has_key(self, option):
        return self.__contains__(option)

    def chk_option(self, option):
        if not option in self.opt_dict:
            raise error_handler.MsgError("Option %s not known!"%option)

    def check_at_lists(self, at_lists, prt_lvl=0):
        """
        Check if an at_lists definition of molecular fragments is useful.
        """
        num_lists = len(at_lists)

        lens = []
        sum_list = []
        for at_list in at_lists:
            sum_list+=at_list
            lens.append(len(at_list))

        numen = len(sum_list)
        maxen = max(sum_list)

        if prt_lvl >= 1:
            print('\nChecking whether the at_lists definition is valid ...')
            if prt_lvl >= 2:
                print('at_lists=', at_lists)
            print('  %i lists with individual numbers of entries:'%(num_lists))
            print(lens)

            print('  %i total entries, with maximal value %i'%(numen,maxen))

        for i in range(1,maxen+1):
            ci = sum_list.count(i)
            if ci!=1:
                print(' WARNING: value %i present %i times in at_lists!'%(i,ci))

    def copy(self, coptions):
        """
        Copy information from a different options instance.
        """
        for key, val in coptions.opt_dict.items():
            self[key] = val

class read_options(options):
    """
    General class for handling input options read from file.
    """
    def __init__(self, ifile, check_init=True):
        options.__init__(self, ifile)

        self.set_defaults()
        self.init = self.read_ifile()

        if check_init: self.check_init()

        self.post_process()

    def check_init(self):
        """
        Check if the instance was properly initialized (the file was read).
        """
        if self.init > 0:
            print("\n ERROR: Input file %s not found!"%self.ifile)
            print("  Please create this file using theoinp")
            exit(0)

    def set_defaults(self):
        """
        Set defaults for the options.
        All possible options should appear here.
        -> inherit for specific implementations
        """
        pass

    def read_ifile(self):
        """
        Read the input file self.ifile.
        Key and value are separated by '='.
        Leading and trailing whitespace is removed.
        """
        if self.ifile is None:
            return 0
        try:
            fileh = open(self.ifile, 'r')
        except:
            return 1

        for line in fileh:
            # take out possible comments
            if '#' in line: continue

            words = line.strip().split('=')
            if len(line.strip()) == 0: continue

            if len(words) != 2:
                print(" ERROR: in file %s\n   line cannot be parsed:"%self.ifile)
                print(len(line))
                print(line)
                exit(6)

            key = words[0].strip()

            if words[1] == '':
                raise error_handler.MsgError('Please specify a value for "%s=" in %s!'%(key, self.ifile))

            val = eval(words[1])

            # every possible option has to be initiliazed in set_defaults to avoid confusion
            if not key in self.opt_dict:
                raise error_handler.MsgError('Unknown option in %s: %s'%(self.ifile, key))

            self.opt_dict[key] = val

        return 0

    def get_def(self, option, default):
        self.chk_option(option)

        if self.opt_dict[option] == None:
            return default
        else:
            return self.opt_dict[option]

    def post_process(self):
        pass

class write_options(options):
    """
    General class for writing options to an input file.
    """
    def __init__(self, ifile):
        options.__init__(self, ifile)

        self.ostr = ''

    def read_str(self, title, key, *args, **kwargs):
        """
        Read a string from input.
        """
        titlek = "%s (%s):"%(title, key)

        val = self.ret_str(titlek, *args, **kwargs)

        self.write_option(key, val)

    def ret_str(self, title, default='', autocomp=False):
        # readline, which is used for auto completion,
        #   creates weird output of the form [?1034h
        #   it should only be imported here
        import readline

        print()
        print(title)

        acstr = ' (autocomplete enabled)' if autocomp else ''
        inpstr = 'Choice%s: '%acstr
        if not default=='': inpstr += '[%s] '%default

        if autocomp:
            readline.set_completer_delims(' \t\n;')
            readline.parse_and_bind("tab: complete")    # activate autocomplete
        val = user_input(inpstr)
        readline.parse_and_bind("tab: ")            # deactivate autocomplete

        if val=='': val = default

        return val

    def read_float(self, title, key, default=1.111):
        """
        Read a float from input.
        """
        titlek = "%s (%s):"%(title, key)

        val = self.ret_float(titlek, default)

        self.write_option(key, val)

    def ret_float(self, title, default=1.111):
        print()
        print(title)

        inpstr = 'Choice: '
        if not default==1.111: inpstr += '[%f] '%default

        sval = user_input(inpstr)
        if sval=='':
            val = default
        else:
            val = float(sval)

        return val

    def read_int(self, title, key, idef=-1):
        """
        Read a string from input.
        """
        titlek = "%s (%s):"%(title, key)

        val = self.ret_int(titlek, idef)

        self.write_option(key, val)

    def ret_int(self, title, idef=-1):
        print()
        print(title)

        return self.inp_int(idef)

    def inp_int(self, idef=-1):
        inpstr = 'Choice: '
        if not idef==-1:
            inpstr += '[%i] '%idef

        retval = idef
        while True:
            try:
                retval = int(user_input(inpstr))
            except:
                if retval==-1:
                    print("Please enter an integer number!")
            if retval!=-1: break

        return retval

    def read_yn(self, title, key, default=False):
        """
        Read Boolean from input.
        """
        titlek = "%s (%s):"%(title, key)

        val = self.ret_yn(titlek, default)

        self.write_option(key, val)

        return val

    def ret_yn(self, question, default=False):
        """
        Ask a yes/no question and return True or False.
        """
        print()
        print(question)

        inpstr = 'Choice (y/n): '
        if default:
            inpstr += '[y] '
        else:
            inpstr += '[n] '

        answer = user_input(inpstr)

        if default:
            return not 'n' in answer.lower()
        else:
            return 'y' in answer.lower()

    def choose_list(self, title, key, opt_expl, default=''):
        """
        Choose an option from a list containing options and explanations.
        """
        titlek = "%s (%s):"%(title, key)

        expl = ["%10s - %s"%(opt, expl) for opt, expl in opt_expl]

        idef = -1
        for ioe, oe in enumerate(opt_expl):
            if oe[0] == default:
                idef = ioe + 1
        ichoice = self.ret_choose_list(titlek, expl, idef)

        val = opt_expl[ichoice-1][0]

        self.write_option(key, val)

    def ret_choose_list(self, title, expl, idef=-1):
        """
        Choose an option from a list containing explanations and return the answer.
        """
        print()
        print(title)

        self.print_list(expl)

        return self.inp_int(idef)

    def print_list(self, plist):
        """
        Print an indexed list to screen.
        """
        iopt = 0
        for p in plist:
            iopt += 1
            print("  [%2i] %s"%(iopt, p))

    def write_list(self, key, wlist, lformat="%i"):
        # write_option can be called directly
        self.write_option(key, wlist)

    def write_option(self, key, val):
        self[key] = val

        if type(val) is str:
            self.ostr += "%s='%s'\n"%(key, str(val))
        else:
            self.ostr += "%s=%s\n"%(key, str(val))

    def flush(self, lvprt=0, choose_file=False):
        if choose_file:
            act_ifile = self.ret_str('Name of input file', self.ifile)
        else:
            act_ifile = self.ifile

        fileh = open(act_ifile, 'w')
        fileh.write(self.ostr)
        fileh.close()
        if lvprt==1:
            print('Finished: File %s written.'%act_ifile)
