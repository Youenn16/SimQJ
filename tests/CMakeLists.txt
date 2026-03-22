add_executable(unit_tests
    test_math.cpp
    test_model.cpp
)

target_link_libraries(unit_tests
    PRIVATE
        core
        Catch2::Catch2WithMain
)

target_include_directories(unit_tests
    PRIVATE
        ../include
)
 
add_test(NAME AllTests COMMAND unit_tests)
catch_discover_tests(unit_tests)

# catch_discover_tests(unit_tests
#     TEST_PREFIX "Catch2_"
#     PROPERTIES LABELS "Catch2"
# )

# Then add manual tag groups with different labels
add_test(NAME MathTagGroup COMMAND unit_tests "[math]")
set_tests_properties(MathTagGroup PROPERTIES LABELS "MathTagGroup")

# Manually add tag-based groups
# add_test(NAME MathTests COMMAND unit_tests "[math]")
# add_test(NAME StringTests COMMAND unit_tests "[model]")


# # Set properties
# set_tests_properties(MathTests PROPERTIES LABELS "[math]")
# set_tests_properties(StringTests PROPERTIES LABELS "[model]")

