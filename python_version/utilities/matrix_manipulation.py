def value_matrix(rows, cols, value):
#    """
#    Creates a matrix filled with the same value.
#        :param rows: the number of rows the matrix should have
#        :param cols: the number of columns the matrix should have
#        :param value: the value to fill the whole matrix

#        :return: list of lists that form the matrix
#    """
    M = []
    while len(M) < rows:
        M.append([])
        while len(M[-1]) < cols:
            M[-1].append(value)

    return M
