## Contributing code

We welcome potential contributions from anyone interested in the project. Any code that you contribute will be 
licensed under the BSD 3-clause license adopted by kmers.

The best way to contribute depends on your expected frequency of contribution.  If believe that you will have 
sustained contributions, or will contribute on multiple occasions over a period of time, then please open an issue 
to request to be added as a contributor to the repository.  If you anticipate making a one time contribution, then 
you can do so directly via a pull request from a fork of the repo. Please make all PRs to the _develop_ branch
of the repository.  PRs made to the _master_ branch may be rejected if they cannot be cleanly rebased
on _develop_.  Before you make a PR, please check that:

 * you've run `cargo fmt` on the relevant code.
 * any non-obvious code is documented (we don't yet have formal documentation guidelines, so use common sense)
 * you've run `cargo clippy` on the relevant code and any issues are either resolved or the PR describes why they were ignored.
