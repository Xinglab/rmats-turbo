from libcpp cimport bool as cbool
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdint cimport *


cdef extern from 'api/BamConstants.h' namespace 'BamTools::Constants' nogil:
    # BAM tag types & sizes
    const char BAM_TAG_TYPE_ASCII
    const char BAM_TAG_TYPE_INT8
    const char BAM_TAG_TYPE_UINT8
    const char BAM_TAG_TYPE_INT16
    const char BAM_TAG_TYPE_UINT16
    const char BAM_TAG_TYPE_INT32
    const char BAM_TAG_TYPE_UINT32
    const char BAM_TAG_TYPE_FLOAT
    const char BAM_TAG_TYPE_STRING
    const char BAM_TAG_TYPE_HEX
    const char BAM_TAG_TYPE_ARRAY


cdef extern from 'api/BamAux.h' namespace 'BamTools' nogil:
    cdef cppclass CigarOp:
        char Type
        uint32_t Length
    
    cdef cppclass RefData:
        string RefName
        int32_t RefLength

    ctypedef vector[RefData] RefVector


cdef extern from 'api/SamHeader.h' namespace 'BamTools' nogil:
    cdef cppclass SamHeader:
        SamHeader()
        SamHeader(const string& headerText)
        SamHeader(const SamHeader& other)

        void Clear()
        string GetErrorString()
        cbool HasError()
        cbool IsValid(cbool verbose = false)
        void SetHeaderText(const string& headerText)
        string ToString()

        cbool HasVersion()
        cbool HasSortOrder()
        cbool HasGroupOrder()
        cbool HasSequences()
        cbool HasReadGroups()
        cbool HasPrograms()
        cbool HasComments()


cdef extern from 'api/BamAlignment.h' namespace 'BamTools' nogil:
    cdef cppclass BamAlignment:
        BamAlignment()
        BamAlignment(const BamAlignment& other)

        cbool IsDuplicate()
        cbool IsFailedQC()
        cbool IsFirstMate()
        cbool IsMapped()
        cbool IsMateMapped()
        cbool IsMateReverseStrand()
        cbool IsPaired()
        cbool IsPrimaryAlignment()
        cbool IsProperPair()
        cbool IsReverseStrand()
        cbool IsSecondMate()

        cbool GetTag[T](const string& tag, T& destination)
        # cbool GetTag[T](const string& tag, vector[T]& destination)

        vector[string] GetTagNames()

        cbool GetTagType(const string& tag, char& type)

        cbool GetArrayTagType(const string& tag, char& type)

        cbool HasTag(const string& tag)

        cbool BuildCharData()

        int GetEndPosition(cbool usePadded = false, cbool closedInterval = false)

        string GetErrorString()

        string Name
        int32_t     Length
        string QueryBases
        string AlignedBases
        string Qualities
        string TagData
        int32_t     RefID
        int32_t     Position
        uint16_t    Bin
        uint16_t    MapQuality
        uint32_t    AlignmentFlag
        vector[CigarOp] CigarData
        int32_t     MateRefID
        int32_t     MatePosition
        int32_t     InsertSize
        string Filename


cdef extern from 'api/BamReader.h' namespace 'BamTools' nogil:
    cdef cppclass BamReader:
        BamReader()

        cbool Close()
        const string GetFilename()
        cbool IsOpen()
        cbool Jump(int refID, int position = 0)
        cbool Open(const string& filename)
        cbool Rewind()

        cbool GetNextAlignment(BamAlignment& alignment)
        cbool GetNextAlignmentCore(BamAlignment& alignment)

        const SamHeader& GetConstSamHeader()
        SamHeader GetHeader()
        string GetHeaderText()

        int GetReferenceCount()
        const RefVector& GetReferenceData()
        int GetReferenceID(const string& refName)

        cbool HasIndex()
        cbool OpenIndex(const string& indexFilename)

        string GetErrorString()
