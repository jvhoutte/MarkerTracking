cmake_minimum_required (VERSION 2.8.12)
project (ImplantPlanning)

find_package (Qt5Widgets)

if(NOT VTK_BINARY_DIR)
	find_package(VTK REQUIRED)
	include(${VTK_USE_FILE})
endif()

set (MyProjectLib_src 
    ${PROJECT_SOURCE_DIR}/mainwindow.cpp 
    ${PROJECT_SOURCE_DIR}/vtkMotionFilter.cpp 
    ${PROJECT_SOURCE_DIR}/vtkAnnotateMarker.cpp 
    ${PROJECT_SOURCE_DIR}/GeometricROI.cpp
    ${PROJECT_SOURCE_DIR}/vtkComputeDistanceField.cpp
    ${PROJECT_SOURCE_DIR}/KneeAPI.cpp
    )
set (MyProjectLib_hdr 
    ${PROJECT_SOURCE_DIR}/mainwindow.h 
    ${PROJECT_SOURCE_DIR}/vtkMotionFilter.h 
    ${PROJECT_SOURCE_DIR}/vtkAnnotateMarker.h
    ${PROJECT_SOURCE_DIR}/GeometricROI.h
    ${PROJECT_SOURCE_DIR}/vtkComputeDistanceField.h
    ${PROJECT_SOURCE_DIR}/KneeAPI.h
)

set (MyProjectLib_ui  ${PROJECT_SOURCE_DIR}/mainwindow.ui)
set (MyProjectBin_src ${PROJECT_SOURCE_DIR}/main.cpp)

qt5_wrap_cpp(MyProjectLib_hdr_moc ${MyProjectLib_hdr})
qt5_wrap_ui (MyProjectLib_ui_moc  ${MyProjectLib_ui})

include_directories (${PROJECT_SOURCE_DIR})
include_directories (${PROJECT_BINARY_DIR})

add_executable(ImplantPlanning
    main.cpp
    ${MyProjectLib_src}
    ${MyProjectLib_hdr_moc}
    ${MyProjectLib_ui_moc}
    ${MyProjectLib_src}
)

target_link_libraries (ImplantPlanning ${VTK_LIBRARIES} Qt5::Widgets )


