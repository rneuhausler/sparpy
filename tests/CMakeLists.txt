set(new_python_path ${CMAKE_BINARY_DIR}:$ENV{PYTHONPATH})
add_custom_target(run_tests 
    PYTHONPATH=${new_python_path} py.test ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Running python tests using pytest..."
    VERBATIM
    )

