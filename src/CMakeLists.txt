include_directories(commons)
add_subdirectory(conterminatorutils)
add_subdirectory(commons)
add_subdirectory(version)
add_subdirectory(workflow)

add_executable(conterminator
        ${commons_source_files}
        ${conterminatorutils_source_files}
        ${workflow_source_files}
        commons/LocalCommandDeclarations.h
        conterminator.cpp
        )

mmseqs_setup_derived_target(conterminator)
target_link_libraries(conterminator version)

add_dependencies(conterminator local-generated)

install(TARGETS conterminator DESTINATION bin)
