Feature: A user can upload a file and save compounds

    Scenario: A flow file can be uploaded for use in compounds
        Given I start the qcluster
        Given I create a project as before
        When I refresh the user object
        Then I can list the projects on the system
        Then the upload URL from the first project in the list points to the right place
        When I upload eight_member_ring.xlsx to flowfiles
        Then The flow file response contains the identifier I gave the file

    Scenario Outline: Data in a file can be previewed based on 20+ examples
        Given I start the qcluster
        Given I create a project as before
        When I refresh the user object
        Then I can list the projects on the system
        Then the upload URL from the first project in the list points to the right place
        When I upload <filename> to flowfiles
        Then The flow file response contains the identifier I gave the file
        When I validate the compounds file
        Then the response from post validate files is accepted

        Examples: CDXML
        |filename|
        |test_cdxml_1.cdxml|


        Examples: SDF
        |filename|
        |test_sdf_1.sdf|
        |test_sdf_2.sdf|
        |test_sdf_3.sdf|
        |test_sdf_4.sdf|
        |test_sdf_5.sdf|
        |test_sdf_6.sdf|
        |test_sdf_7.sdf|
        |test_sdf_8.sdf|
        |file_1.sdf|

        Examples: XLSX
        |filename|
        |test_excel_1.xlsx|
        |test_excel_2.xlsx|
        |test_excel_3.xlsx|
        |test_excel_4.xlsx|
        |test_excel_5.xlsx|
        |test_excel_6.xlsx|
        |eight_member_ring.xlsx|



