
Feature: Upload contains stereochemistry
'''
    I wish to submit a compound or compounds which contain stereochemistry and for that stereochemistry to be retained.

    Scenario: User submits a substance containing one stereocentre to proja. 
        Given a User
        and my user is member of a group
        and a valid project exists proja
        and I have valid stereochem compound with one stereocentre
        #and I have a valid molfile
        and I automatically have editor permissions as creator
        When I log in 
        then I validate my cbh_compound_batch to proja
        and I create my cbh_compound_batch to proja
        and retain its stereochemistry

    Scenario: User submits a substance containing multiple stereocentres to proja. 
        Given a User
        and my user is member of a group
        and a valid project exists proja
        and I have valid stereochem compound with multiple stereocentres
        #and I have a valid molfile
        and I automatically have editor permissions as creator
        When I log in 
        then I validate my cbh_compound_batch to proja
        and I create my cbh_compound_batch to proja
        and retain its stereochemistry

    Scenario: User submits a Chemdraw file containing multiple compounds which contain stereochemistry to proja. 
        Given a User
        and my user is member of a group
        and a valid project exists proja
        and I have valid stereochem compounds within a ChemDraw file
        and I automatically have editor permissions as creator
        When I log in 
        then I validate my stereochem compounds to proja
        and I create my stereochem compounds to proja
        and retain all their stereochemistry

	Scenario: User submits a SD file containing multiple compounds which contain stereochemistry to proja. 
	    Given a User
        and my user is member of a group
        and a valid project exists proja
        and I have valid stereochem compounds within a SD file
        and I automatically have editor permissions as creator
        When I log in 
        then I validate my stereochem compounds to proja
        and I create my stereochem compounds to proja
        and retain all their stereochemistry

