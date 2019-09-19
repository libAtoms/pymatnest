import os, psutil, gc

active=False
last_rss = 0
def check_memory(label, do_garbage=False):
    global last_rss, active

    if not active:
        return

    gc.collect()

    process = psutil.Process(os.getpid())
    m = process.memory_info()
    print label, "MEMORY", m.rss-last_rss, m.rss, m.vms  # in bytes 
    if do_garbage:
        print label, "GARBAGE", gc.garbage
    last_rss = m.rss
