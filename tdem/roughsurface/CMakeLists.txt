ADD_EXECUTABLE        (roughsurface   roughsurface.cpp)
TARGET_LINK_LIBRARIES (roughsurface   ${LIBS})
SET_TARGET_PROPERTIES (roughsurface   PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
