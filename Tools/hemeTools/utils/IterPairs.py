def IterPairs(iterable):
    """Return adjacent pairs of items from iterable. I.e.:
    IterPairs(lst) -> (lst[0], lst[1]), (lst[1], lst[2]), (lst[2], lst[3]), ..., (lst[i], lst[i+1])
    """
    first = True
    for item in iterable:
        if first:
            new = item
            first = False
            continue
        else:
            old = new
            new = item
            yield old, new
        continue
    return
