'use strict';

describe('Service: CBHCompoundBatch', function () {

  // load the service's module
  beforeEach(module('ngChemApp'));

  // instantiate service
  var CBHCompoundBatch;
  beforeEach(inject(function (_CBHCompoundBatch_) {
    CBHCompoundBatch = _CBHCompoundBatch_;
  }));

  it('should do something', function () {
    expect(!!CBHCompoundBatch).toBe(true);
  });

});
