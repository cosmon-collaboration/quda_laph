add_executable(TestRephase test_rephase.cc)
target_include_directories(TestRephase PUBLIC include)
target_link_libraries(TestRephase PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)

add_executable(TestMomenta test_Momenta.cc)
target_include_directories(TestMomenta PUBLIC include)
target_link_libraries(TestMomenta PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)

add_executable(TestPeramSrc test_peram_src.cc)
target_include_directories(TestPeramSrc PUBLIC include)
target_link_libraries(TestPeramSrc PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)

add_executable(TestNoise test_noise.cc)
target_include_directories(TestNoise PUBLIC include)
target_link_libraries(TestNoise PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)

add_executable(TestColorCross test_ColorCross.cc)
target_include_directories(TestColorCross PUBLIC include)
target_link_libraries(TestColorCross PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)

add_executable(TestColorContract test_ColorContract.cc)
target_include_directories(TestColorContract PUBLIC include)
target_link_libraries(TestColorContract PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)

add_executable(TestBaryonKernel test_laphBaryonKernel.cc)
target_include_directories(TestBaryonKernel PUBLIC include)
target_link_libraries(TestBaryonKernel PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)

add_executable(TestBaryonTripletA test_laphBaryonTripletA.cc)
target_include_directories(TestBaryonTripletA PUBLIC include)
target_link_libraries(TestBaryonTripletA PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)

add_executable(TestBaryonTripletB test_laphBaryonTripletB.cc)
target_include_directories(TestBaryonTripletB PUBLIC include)
target_link_libraries(TestBaryonTripletB PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)

add_executable(TestMesonDoublet test_laphMesonDoublet.cc)
target_include_directories(TestMesonDoublet PUBLIC include)
target_link_libraries(TestMesonDoublet PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)

add_executable(TestBatchGEMM test_BatchGEMM.cc)
target_include_directories(TestBatchGEMM PUBLIC include)
target_link_libraries(TestBatchGEMM PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)

add_executable(TestMesonKernel test_laphMesonKernel.cc)
target_include_directories(TestMesonKernel PUBLIC include)
target_link_libraries(TestMesonKernel PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)

add_executable(TestQuartet test_laphQuartet.cc)
target_include_directories(TestQuartet PUBLIC include)
target_link_libraries(TestQuartet PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)

add_executable(TestQuartetKernel test_laphQuartetKernel.cc)
target_include_directories(TestQuartetKernel PUBLIC include)
target_link_libraries(TestQuartetKernel PUBLIC ${LiveLaphLibs} gputoolkitlib quda QudaLaph)
