'use strict';

describe('Service: chemFileValidation', function () {

  // load the service's module
  beforeEach(module('ngChemApp'));

  // instantiate service
  var chemFileValidation;
  beforeEach(inject(function (_chemFileValidation_) {
    chemFileValidation = _chemFileValidation_;
  }));

  it('should do something', function () {
    expect(!!chemFileValidation).toBe(true);
  });

});
