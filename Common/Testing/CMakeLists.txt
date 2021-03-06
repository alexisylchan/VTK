ADD_SUBDIRECTORY(Cxx)

IF (VTK_WRAP_TCL)
  ADD_SUBDIRECTORY(Tcl)
ENDIF (VTK_WRAP_TCL)

IF (VTK_WRAP_PYTHON)
  ADD_SUBDIRECTORY(Python)
ENDIF (VTK_WRAP_PYTHON)

IF(PYTHON_EXECUTABLE)
  ADD_TEST(HeaderTesting-Common ${PYTHON_EXECUTABLE}
    ${VTK_SOURCE_DIR}/Common/Testing/HeaderTesting.py
    "${VTK_SOURCE_DIR}/Common"
    VTK_COMMON_EXPORT
    vtkABI.h
    vtkAbstractIterator.h
    vtkAbstractList.h
    vtkAbstractMap.h
    vtkArray.h
    vtkArrayCoordinateIterator.h
    vtkArrayCoordinates.h
    vtkArrayExtents.h
    vtkArrayExtentsList.h
    vtkArrayInterpolate.h
    vtkArrayIteratorIncludes.h
    vtkArrayIteratorTemplate.h
    vtkArrayMap.h
    vtkArrayMapIterator.h
    vtkArrayPrint.h
    vtkArrayRange.h
    vtkArraySort.h
    vtkArrayWeights.h
    vtkAssemblyPaths.h
    vtkBoundingBox.h
    vtkBreakPoint.h
    vtkByteSwap.h
    vtkCallbackCommand.h
    vtkCommand.h
    vtkCommonInformationKeyManager.h
    vtkContainer.h
    vtkDataArrayCollection.h
    vtkDataArrayTemplate.h
    vtkDebugLeaks.h
    vtkDebugLeaksManager.h
    vtkDenseArray.h
    vtkDynamicLoader.h
    vtkErrorCode.h
    vtkEventForwarderCommand.h
    vtkFloatingPointExceptions.h
    vtkGarbageCollector.h
    vtkGarbageCollectorManager.h
    vtkHashMap.h
    vtkHashMapIterator.h
    vtkIOStream.h
    vtkIOStreamFwd.h
    vtkIdListCollection.h
    vtkImplicitFunctionCollection.h
    vtkIndent.h
    vtkInformation.h
    vtkInformationInternals.h
    vtkJavaAwt.h
    vtkJavaUtil.h
    vtkLargeInteger.h
    vtkLinkedList.h
    vtkLinkedListIterator.h
    vtkNew.h
    vtkOStrStreamWrapper.h
    vtkOStreamWrapper.h
    vtkObject.h
    vtkObjectBase.h
    vtkObjectFactoryCollection.h
    vtkOldStyleCallbackCommand.h
    vtkOnePieceExtentTranslator.h
    vtkOverrideInformationCollection.h
    vtkPlaneCollection.h
    vtkProcessStatistics.h
    vtkPropCollection.h
    vtkPythonCommand.h
    vtkRayCastStructures.h
    vtkRungeKutta2.h 
    vtkSetGet.h
    vtkSmartPointer.h
    vtkSmartPointerBase.h
    vtkSparseArray.h
    vtkStdString.h
    vtkStructuredData.h
    vtkSystemIncludes.h
    vtkTcl.h
    vtkTclUtil.h
    vtkTemplateAliasMacro.h
    vtkTimeStamp.h
    vtkTk.h
    vtkTransformCollection.h
    vtkTryDowncast.h
    vtkType.h
    vtkTypeTemplate.h
    vtkTypeTraits.h
    vtkTypedArray.h
    vtkUnicodeString.h
    vtkVariant.h
    vtkVariantBoostSerialization.h
    vtkVariantCast.h
    vtkVariantCreate.h
    vtkVariantExtract.h
    vtkVariantInlineOperators.h
    vtkVector.h
    vtkVectorIterator.h
    vtkVersion.h
    vtkWeakPointer.h
    vtkWeakPointerBase.h
    vtkWin32Header.h
    vtkWindows.h
    vtkXMLFileOutputWindow.h
    )
ENDIF(PYTHON_EXECUTABLE)

#-----------------------------------------------------------------------------
# Configure the VTK install test.  It is enabled only if the install
# prefix is set to an InstallTest/Root directory under the top of the build
# tree.  Note that the install prefix used for the install test must
# not be a 1-level deep subdirectory under the top of the build tree.
# If it were such then a few executables may be confused into thinking
# they are running from the build tree instead of the install tree.
IF(WIN32 OR APPLE)
  STRING(TOLOWER "${CMAKE_INSTALL_PREFIX}" VTK_TEST_INSTALL_LEFT)
  STRING(TOLOWER "${VTK_BINARY_DIR}/InstallTest/Root" VTK_TEST_INSTALL_RIGHT)
ELSE(WIN32 OR APPLE)
  SET(VTK_TEST_INSTALL_LEFT "${CMAKE_INSTALL_PREFIX}")
  SET(VTK_TEST_INSTALL_RIGHT "${VTK_BINARY_DIR}/InstallTest/Root")
ENDIF(WIN32 OR APPLE)
STRING(COMPARE EQUAL "${VTK_TEST_INSTALL_LEFT}" "${VTK_TEST_INSTALL_RIGHT}"
  VTK_TEST_INSTALL)

# We cannot safely implement the install test unless the configuration
# type is known when the test is run.
IF(NOT VTK_TEST_CONFIG_TYPE_KNOWN)
  SET(VTK_TEST_INSTALL 0)
ENDIF(NOT VTK_TEST_CONFIG_TYPE_KNOWN)

IF(VTK_TEST_INSTALL)
  IF(CMAKE_CONFIGURATION_TYPES)
    # There are multiple configurations.  Make sure the tested
    # configuration is the one that is installed.
    SET(DOLLAR "$")
    SET(VTK_INSTALL_TEST_CONFIG_TYPE -C "${DOLLAR}{CTEST_CONFIGURATION_TYPE}")
  ELSE(CMAKE_CONFIGURATION_TYPES)
    # There is only one configuration.  It will be installed.
    SET(VTK_INSTALL_TEST_CONFIG_TYPE)
  ENDIF(CMAKE_CONFIGURATION_TYPES)

  # Add a test to install VTK through the build system install target.
  ADD_TEST(Install
    ${CMAKE_CTEST_COMMAND}
    ${VTK_INSTALL_TEST_CONFIG_TYPE}
    --build-and-test ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR}
    --build-generator ${CMAKE_GENERATOR}
    --build-project VTK
    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
    --build-noclean
    --build-target install)
ENDIF(VTK_TEST_INSTALL)

# Suppress memory checking of some tests
CONFIGURE_FILE(${VTK_SOURCE_DIR}/Common/Testing/CTestCustom.ctest.in
  ${VTK_BINARY_DIR}/Common/Testing/CTestCustom.ctest @ONLY)
