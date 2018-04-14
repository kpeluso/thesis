# Source:
#   https://stackoverflow.com/questions/492519/timeout-on-a-function-call

def timeout_inny(secs, fun, *params):
    '''
    Function 'fun' with arguments 'params' ~may~ run for 'secs' seconds.
    '''
    import signal

    # Register an handler for the timeout
    def handler(signum, frame):
        raise Exception("FROM 'timeout(): fun has runout of 'secs' time!")

    signal.signal(signal.SIGALRM, handler) # Register the signal function handler
    signal.alarm(secs) # Define a timeout for function

    # Run function until timeout
    try:
        fun(*params)
        return 1
    except Exception, exc:
        print exc
        return 0


def timeout(secs, fun, *params):
    '''
    Function 'fun' with arguments 'params' ~may~ run for 'secs' seconds.
    '''
    import signal
    if timeout_inny(secs, fun, *params) == 1:
        signal.alarm(0) # end signal thread!
        return 1
    else:
        return 0

# FOR TESTING:
# from time import sleep
# from TIMEOUT_clean import timeout
# def t(n):
#     sleep(n)
#     return 0
# timeout(5,t,2)
# timeout(2,t,5)

