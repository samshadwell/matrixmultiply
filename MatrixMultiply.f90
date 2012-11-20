!Very rudimentary Fox Matrix Multiplication program
!@author Samuel Shadwell
!
!KNOWN ISSUES: Only works with square matrices
!   Requires more than 1 MPI process


PROGRAM MatrixMultiply
IMPLICIT NONE
INCLUDE 'mpif.h'

!Variable declarations
INTEGER                                 :: ierr, comm, commSize, rank           !MPI variables
REAL, ALLOCATABLE, DIMENSION(:,:)       :: arayA, arayB, arayC                  !Matrices
INTEGER                                 :: arayLen
REAL                                    :: nodeAVal, nodeBVal

!Initialize MPI Environment
CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE (comm, commSize, ierr)
CALL MPI_COMM_RANK (comm, rank, ierr)

!Sets and distributes array values to various nodes
CALL setUpArrays(ierr, comm, commSize, rank, arayA, arayB, arayC, nodeAVal, nodeBVal, arayLen)

!Parallel Code Ends
CALL MPI_FINALIZE (ierr)

!Final Serial Code


END PROGRAM MatrixMultiply



!============== SUBROUTINE SetUpArrays ===============
!   Puts entire array on head node (0), and distributes single elements to 
!   all other nodes. Also sets size of arayA, B, and C to 1 for all
!   MPI processes other than the head node.
SUBROUTINE setUpArrays(ierr, comm, commSize, rank, arayA, arayB, arayC, nodeAVal, nodeBVal, arayLen)
IMPLICIT NONE
INTEGER                                 :: ierr, comm, commSize, rank, arayLen, mpiStatus
REAL, ALLOCATABLE, DIMENSION(:,:)       :: arayA, arayB, arayC
REAL                                    :: nodeAVal, nodeBVal
INTEGER                                 :: i, j, n, tag


!Currently sets values for arrays, I need to figure out how to read them in

arayLen = 3                 ! Array dimensions are 3x3

IF (rank .EQ. 0) THEN       ! Begin head node

    ALLOCATE(arayA(arayLen, arayLen), arayB(arayLen, arayLen), arayC(arayLen,arayLen))
    
    !Populates arrays A and B
    !TODO: replace this with an actual parser
    DO i = 1, arayLen
        DO j = 1, arayLen
            arayA(i,j) = i+j
            arayB(i,j) = i-j
            arayC(i,j) = 0  !Initialize all C elements as 0
        ENDDO
    ENDDO

    !Sends all nodes their values
    n = arayLen**2          !Number of data values which must be sent
    DO i = 1, n
        CALL MPI_SEND(arayA(MODULO(i, arayLen), i/arayLen), 1, MPI_REAL, MODULO(i, commSize-1), 1, comm, ierr)
        CALL MPI_SEND(arayB(MODULO(i, arayLen), i/arayLen), 1, MPI_REAL, MODULO(i, commSize-1), 2, comm, ierr)
    ENDDO

ENDIF                       ! End head node


IF (rank .NE. 0) THEN       ! Nodes other than head node

    ALLOCATE(arayA(1,1), arayB(1,1), arayC(1,1))
    !Sets C array to 0
    arayC(1,1) = 0.0
    !Gets value for arrays A and B sent from head node
    mpiStatus = status(MPI_STATUS_SIZE)

    CALL MPI_RECV(nodeAVal, 1, MPI_REAL, 0, 1, comm, mpiStatus, ierr)
    CALL MPI_RECV(nodeBVal, 1, MPI_REAL, 0, 2, comm, mpiStatus, ierr)


ENDIF                       ! End other nodes

END SUBROUTINE setUpArrays
