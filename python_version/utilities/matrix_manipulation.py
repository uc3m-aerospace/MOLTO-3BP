def value_matrix(rows, cols, value):
#    """
#    Creates a matrix filled with the same value.
#        :param rows: the number of rows the matrix should have
#        :param cols: the number of columns the matrix should have
#        :param value: the value to fill the whole matrix

#        :return: list of lists that form the matrix
#    """
    M = []
    if rows > 1:
        while len(M) < rows:
            M.append([])
            while len(M[-1]) < cols:
                M[-1].append(value)
    else:
        while len(M) < cols:
            M.append(value)

    return M

def rand_value_matrix(rows, cols):
#    """
#    Creates a matrix filled with the random value.
#        :param rows: the number of rows the matrix should have
#        :param cols: the number of columns the matrix should have

#        :return: list of lists that form the matrix
#    """
    import random
    M = []
    if rows > 1:
        while len(M) < rows:
            M.append([])
            while len(M[-1]) < cols:
                M[-1].append(random.random())
    else:
        while len(M) < cols:
            M.append(random.random())

    return M

def matrix_multiply(A, B):
#    """
#    Returns the product of the matrix A * B
#        :param A: The first matrix - ORDER MATTERS!
#        :param B: The second matrix

#        :return: The product of the two matrices
#    """
    # Section 1: Ensure A & B dimensions are correct for multiplication
    rowsA = len(A)
    if isinstance(A[0],list):
        colsA = len(A[0])
    else:
        colsA = 1
    rowsB = len(B)
    if isinstance(B[0],list):
        colsB = len(B[0])
    else:
        colsB = 1
    if colsA != rowsB:
        raise ArithmeticError(
            'Number of A columns must equal number of B rows.')

    # Section 2: Store matrix multiplication in a new matrix
    C = value_matrix(rowsA, colsB, 0)
    for i in range(rowsA):
        for j in range(colsB):
            total = 0
            for ii in range(colsA):
                total += A[i][ii] * B[ii][j]
            C[i][j] = total

    return C

def vector_multiply(A, B):
#    """
#    Returns the product of the vectors A * B
#        :param A: The first row vector - ORDER MATTERS!
#        :param B: The second row vector

#        :return: The product of the two vectors (Matrix)
#    """
    # Section 1: Take A & B dimensions
    colsA = len(A)
    colsB = len(B)

    # Section 2: Store matrix multiplication in a new matrix
    C = value_matrix(colsA, colsB, 0)
    for i in range(colsA):
        for j in range(colsB):
            C[i][j] = A[i]*B[j]

    return C

def transpose(M):
#    """
#    Returns a transpose of a matrix.
#        :param M: The matrix to be transposed
#
#        :return: The transpose of the given matrix
#    """
    # Section 1: if a 1D array, convert to a 2D array = matrix
    if not isinstance(M[0],list):
        M = [M]

    # Section 2: Get dimensions
    rows = len(M)
    cols = len(M[0])

    # Section 3: MT is zeros matrix with transposed dimensions
    MT = value_matrix(cols, rows, 0)

    # Section 4: Copy values from M to it's transpose MT
    for i in range(rows):
        for j in range(cols):
            MT[j][i] = M[i][j]

    return MT

def matrix_subtraction(A, B):
#    """
#    Subtracts matrix B from matrix A and returns difference
#        :param A: The first matrix
#        :param B: The second matrix

#        :return: Matrix difference
#    """
    # Section 1: Ensure dimensions are valid for matrix subtraction
    rowsA = len(A)
    if isinstance(A[0],list):
        colsA = len(A[0])
    else:
        colsA = 1
    rowsB = len(B)
    if isinstance(B[0],list):
        colsB = len(B[0])
    else:
        colsB = 1
    if rowsA != rowsB or colsA != colsB:
        raise ArithmeticError('Matrices are NOT the same size.')

    # Section 2: Create a new matrix for the matrix difference
    C = value_matrix(rowsA, colsB, 0)

    # Section 3: Perform element by element subtraction
    if isinstance(A[0],list) and isinstance(B[0],list):
        for i in range(rowsA):
            for j in range(colsB):
                C[i][j] = A[i][j] - B[i][j]
    else:
        for i in range(rowsA):
            C[i] = A[i] - B[i]

    return C
