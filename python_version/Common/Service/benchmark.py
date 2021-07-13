import inspect
import os
import resource
import time
from threading import Thread

from Common.Base.singleton import SingletonMeta
from Common.Config.configurator import Configurator


if Configurator().get('logs.use_colors', False):
    class BColors:
        HEADER = '\033[95m'
        OKBLUE = '\033[94m'
        OKCYAN = '\033[96m'
        OKGREEN = '\033[92m'
        WARNING = '\033[93m'
        FAIL = '\033[91m'
        ENDC = '\033[0m'
        BOLD = '\033[1m'
        UNDERLINE = '\033[4m'
else:
    class BColors:
        HEADER = ' '
        OKBLUE = ' '
        OKCYAN = ' '
        OKGREEN = ' '
        WARNING = ' '
        FAIL = ' '
        ENDC = ' '
        BOLD = ' '
        UNDERLINE = ' '


_ts = time.time()
_debug_trace_string = f"{BColors.BOLD}[%s]{BColors.ENDC} %s ({BColors.OKGREEN}%.3fMb {BColors.ENDC}-{BColors.OKCYAN} %.4fs{BColors.ENDC})"


class Benchmark(metaclass=SingletonMeta):
    """
    Class to trace memory and time for python code lines
    """
    ts_stack = []
    mem_stack = []
    func_stack = []
    key_stack = []

    def reset_time(self):
        self.clear_stacks()
        global _ts
        _ts = time.time()

    def clear_stacks(self):
        """
        Clear all stacks
        """
        self.ts_stack.clear()
        self.mem_stack.clear()
        self.func_stack.clear()
        self.key_stack.clear()

    def _get_memory_usage(self):
        try:
            usage = resource.getrusage(resource.RUSAGE_SELF)
            memory_used = round(float(usage[2]) / (1024.0 * 1024.0), 3)
        except:
            memory_used = 0.0  # non-Linux?

        self.mem_stack.append(memory_used)
        return memory_used

    def _get_timestamp(self):
        global _ts
        _currentTs = time.time()
        timestamp = round(_currentTs - _ts, 4)
        self.ts_stack.append(timestamp)
        return timestamp

    def _get_file_called(self, relative_frame=1):
        total_stack = inspect.stack()
        frame_info = total_stack[relative_frame][0]
        frame = "%s [%d]" % \
               (os.path.basename(frame_info.f_code.co_filename), frame_info.f_lineno)
        self.func_stack.append(frame)
        return frame

    def _get_parsed_stacks(self):
        stack = []
        for i in range(len(self.func_stack)):
            stack.append(_debug_trace_string % (self.key_stack[i], self.func_stack[i], self.mem_stack[i], self.ts_stack[i]))
        return stack

    @staticmethod
    def get_memory_usage():
        """
        Method to have the current memory usage
        """
        return Benchmark()._get_memory_usage()

    @staticmethod
    def get_timestamp():
        """
        Method to gather the elapsed time from the start of execution to now
        """
        return Benchmark()._get_timestamp()

    @staticmethod
    def get_called(relative_frame=1):
        """
        Method to gather the filename and the line executed
        """
        return Benchmark()._get_file_called(relative_frame)

    @staticmethod
    def stat(key="Log", print_output=False):
        """
        Method to add a new trace to the benchmark stack
        """
        trace = _debug_trace_string % \
                (key, Benchmark.get_called(3), Benchmark.get_memory_usage(),
                 Benchmark.get_timestamp())
        Benchmark().key_stack.append(key)
        if print_output:
            print(trace)

    @staticmethod
    def print_stack(clear_stacks=False):
        """
        Method to print the complete stack
        """
        for message in Benchmark()._get_parsed_stacks():
            print(message)

        if clear_stacks:
            Benchmark().clear_stacks()


if __name__ == "__main__":
    """
    Example of use
    """
    Benchmark.stat('Start benchmark', True) # Trace and print
    print('-----------')
    array = []
    for i in range(25):
        time.sleep(0.001)
        Benchmark.stat('First Iteration %d' % i)
    Benchmark.print_stack(True) # Print all stacks and clear the current traces

    print('-----------')
    for i in range(25):
        time.sleep(0.001)
        Benchmark.stat('Second Iteration %d' % i)
    Benchmark.stat('Finished benchmark')  # Trace without print
    Benchmark.print_stack(True)

    print('-----------')

    def test_thread(test_number=1):
        time.sleep(0.1)
        Benchmark.stat("Thread %d" % test_number)

    Benchmark.stat('Start threads')  # Trace in thread mode
    trace1 = Thread(target=test_thread, args=(1,))
    trace2 = Thread(target=test_thread, args=(2,))
    trace3 = Thread(target=test_thread, args=(3,))
    trace1.run()
    trace2.run()
    trace3.run()
    Benchmark.stat('End threads')  # Trace in thread mode
    Benchmark.print_stack(True)
