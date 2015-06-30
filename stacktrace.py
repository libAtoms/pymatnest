# from stack exchange

import code, traceback, signal, sys

def debug(sig, frame):
    """Interrupt running process, and provide a python prompt for
    interactive debugging."""
    d={'_frame':frame}         # Allow access to frame object.
    d.update(frame.f_globals)  # Unless shadowed by global
    d.update(frame.f_locals)

    # i = code.InteractiveConsole(d)
    message  = "Signal received : entering python shell.\nTraceback:\n"
    message += ''.join(traceback.format_stack(frame))
    sys.stderr.write(message+"\n")
    # i.interact(message)
    sys.exit(1)

def listen():
    signal.signal(signal.SIGUSR1, debug)  # Register handler
