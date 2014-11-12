class DelaysFileWrong(Exception):
    pass


class TooManyCandidates(Exception):
    def __init__(self, n):
        self.msg = 'Data set contains too many candidates max = %d !' % n

    def __str__(self):
        return self.msg

