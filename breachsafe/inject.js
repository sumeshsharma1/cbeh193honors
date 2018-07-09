HTMLFormElement.prototype._nativeSubmit = HTMLFormElement.prototype.submit;
HTMLFormElement.prototype.submit = function () {
  var submitEvent = document.createEvent("HTMLEvents");
  submitEvent.initEvent("submit", true, true);
  if (this.dispatchEvent(submitEvent)) {
    this._nativeSubmit.apply(this, arguments);
  }
};
