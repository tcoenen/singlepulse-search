def check_delays_option(options, args, p):
    '''
    Check --delays commandline option.
    '''
    if options.delays is not None:
        delays_file = os.path.realpath(options.delays)
        if not os.path.exists(delays_file):
            print 'Delays file %s does not exist!' % delays_file
            p.print_help()
            sys.exit(1)
    else:
        delays_file = ''

    return delays_file

