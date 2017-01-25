
# [RMS] TODO: need to make this work for xcode?

# set PCH for VS project
function(SET_TARGET_PRECOMPILED_HEADER Target PrecompiledHeader PrecompiledSource)
  if(MSVC)
     SET_TARGET_PROPERTIES(${Target} PROPERTIES COMPILE_FLAGS "/Yu${PrecompiledHeader}")
     set_source_files_properties(${PrecompiledSource} PROPERTIES COMPILE_FLAGS "/Yc${PrecompiledHeader}")
  endif(MSVC)
endfunction(SET_TARGET_PRECOMPILED_HEADER)