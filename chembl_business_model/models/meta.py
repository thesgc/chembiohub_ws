__author__ = 'mnowotka'

from django.db import models
from chembl_core_db.db.models.abstractModel import ChemblAppAbstractModel

#-----------------------------------------------------------------------------------------------------------------------

class InchiErrors(ChemblAppAbstractModel):
    #api_exclude = []

    ERROR_TYPE_CHOICES = (
        ('1.02', 'Inchi Software version 1.02'),
        ('1.03', 'Inchi Software version 1.03'),
        ('1.04', 'Inchi Software version 1.04'),
        ('1.04 RD', 'RDkit using Inchi ver. 1.04'),
        ('1.03 RD', 'RDkit using Inchi ver. 10,3'),
        ('mol read 1.4', 'Read molfile by RDkit using inchi 1.04'),
        ('mol read 1.3', 'Read molfile by RDkit using inchi 1.03'),
        ('open babel 2.3.2', 'Open Babel software version 2.3.2'),
        ('openbabel 2.3.2 runtime', 'Open Babel 2.3.2 runtime error'),
        ('indigo 1.1.5.0 linux32', 'Indigo software, version 1.1.5.0'),
        ('indigo 1.1.5.0 runtime', 'Indigo software, version 1.1.5.0 runtime error'),
        )

    error_type = models.CharField(max_length=30, choices=ERROR_TYPE_CHOICES)
    structure = models.ForeignKey('chembl_core_model.CompoundStructures')

    class Meta(ChemblAppAbstractModel.Meta):
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class ImageErrors(ChemblAppAbstractModel):
    #api_exclude = []

    ERROR_TYPE_CHOICES = (
        ('reg', 'Regular (big) image'),
        ('tomb', 'Thumbnail'),
        )

    error_type = models.CharField(max_length=10, choices=ERROR_TYPE_CHOICES)
    image = models.ForeignKey('chembl_core_model.CompoundImages')

    class Meta(ChemblAppAbstractModel.Meta):
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class DjangoCheatSheet(ChemblAppAbstractModel):
    bigIntegerField = models.BigIntegerField()
    booleanField = models.BooleanField()
    charField = models.CharField(max_length=10)
    commaSeparatedIntegerField = models.CommaSeparatedIntegerField(max_length=20)
    dateField = models.DateField()
    dateTimeField = models.DateTimeField()
    decimalField = models.DecimalField(max_digits=9, decimal_places=3)
    emailField = models.EmailField()
    filePathField = models.FilePathField()
    floatField = models.FloatField()
    integerField = models.IntegerField()
    ipAddressField = models.GenericIPAddressField()
    genericIPAddressField = models.GenericIPAddressField()
    nullBooleanField = models.NullBooleanField()
    positiveIntegerField = models.PositiveIntegerField()
    positiveSmallIntegerField = models.PositiveSmallIntegerField()
    slugField = models.SlugField()
    smallIntegerField = models.SmallIntegerField()
    textField = models.TextField()
    timeField = models.TimeField()
    urlField = models.URLField()

    class Meta(ChemblAppAbstractModel.Meta):
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------

class SDF(ChemblAppAbstractModel):
    originalSDF = models.FileField(upload_to='sdfs')
    originalHash = models.CharField(unique=True, max_length=32)
    cleanSDF = models.FileField(upload_to='sdfs')
    cleanHash = models.CharField(primary_key=True, max_length=32)

    class Meta(ChemblAppAbstractModel.Meta):
        app_label = 'chembl_business_model'

#-----------------------------------------------------------------------------------------------------------------------