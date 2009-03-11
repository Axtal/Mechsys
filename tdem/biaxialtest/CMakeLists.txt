ADD_EXECUTABLE        (biaxialtest   biaxialtest.cpp)
TARGET_LINK_LIBRARIES (biaxialtest   ${LIBS})
SET_TARGET_PROPERTIES (biaxialtest   PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
