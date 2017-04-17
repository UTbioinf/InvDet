#!/usr/bin/env python

import os
import time

class StagingException(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

class ExecStaging(object):
    """Staging small tasks of a large project, so that the whole project can be resumed in the
    middle. The tasks can have stack-like levels
    """
    def __init__(self, directory = ""):
        self._fname = os.path.join(directory, "staging")
        self._write = None
        self._embed_cnt = 0
    
    def resume(self):
        """Resume the task:
        Return values: a dictionary {
                "event_stack": <the unfinished event stack>,
                "cur_unfinished_event": <the last unfinished event in the deepest stack level>
                "last_finished_event": <the previous finished event in the deepest stack levle>
                }
        To resume the task, do the following:
            Follow the "event_stack" to the deepest unfinished stack level
            if cur_unfinished_event != None:
                start from there
            else:
                start from the task after `last_finished_event`
        """
        if not os.path.isfile(self._fname):
            return {"event_stack": [], 
                    "cur_unfinished_event": None, 
                    "last_finished_event": None}
        with open(self._fname, "r") as fin:
            unfinished_stack = []
            last_finished_task = None
            cur_unfinished_task = None
            while True:
                line = fin.readline()
                if not line: break
                tokens = line.split(';', 1)
                if tokens[0][-1] == 's':
                    self._embed_cnt += 1
                    cur_unfinished_task = tokens[1].strip()
                    unfinished_stack.append( cur_unfinished_task )
                else:
                    self._embed_cnt -= 1
                    if self._embed_cnt < 0:
                        raise StagingException("Event not match! See '{}' for details".format(self._fname))
                    last_finished_task = unfinished_stack.pop()
                    cur_unfinished_task = None
                    
            return {"event_stack": unfinished_stack, 
                    "cur_unfinished_event": cur_unfinished_task, 
                    "last_finished_event": last_finished_task}

            

    def begin(self):
        """Start staging from the last aborted task, if any"""
        self._write = open(self._fname, "a")

    def startover(self):
        """Start staging from beginning"""
        self._write = open(self._fname, "w")

    def end(self):
        """Stop staging"""
        self._write.close()

    def abort(self):
        """Abort staging, same as end()"""
        self._write.close()

    def push_event(self, event):
        cur_time = time.time()
        self._write.write("{}{},{},s; {}\n".format(
                '..'*self._embed_cnt,
                cur_time, 
                time.ctime( cur_time ), 
                event))
        self._embed_cnt += 1

    def pop_event(self):
        cur_time = time.time()
        self._embed_cnt -= 1
        if self._embed_cnt < 0:
            raise StagingException("Event not match! See '{}' for details".format(self._fname))
        self._write.write("{}{},{},e;\n".format(
                '..'*self._embed_cnt,
                cur_time, 
                time.ctime( cur_time )
                ))

def test():
    import pprint
    pp = pprint.PrettyPrinter(indent = 4)
    staging = ExecStaging()
    staging.startover()
    staging.push_event("1")
    staging.push_event("1.1")
    staging.pop_event()
    staging.push_event("1.2")
    staging.push_event("1.2.1")
    staging.pop_event()
    staging.push_event("1.2.2")
    staging.pop_event()
    staging.abort()
    pp.pprint( staging.resume() )

    staging.begin()
    staging.push_event("1.2.3")
    staging.abort()
    pp.pprint( staging.resume() )
    

if __name__ == "__main__":
    test()
