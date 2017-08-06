#pragma once

#include <memory>
#include <string>
#include "utility/proc_data.hpp"
#include "utility/thread_squad.hpp"

namespace pipeline {

class AbstractPipelineHead;

class AbstractPipelineNode
{
  /** Generally accepted convention is that the 'friend' designation should be
   * used sparingly, but in consideration of the special relationship between
   * the pipeline head and nodes, its use here is tolerated to ease design
   * and avoid sacrificing access appropriate access limitations.
   * Specifically, AbstractPipelineNode needs to invoke
   * AbstractPipelineHead::SubmitThreadWorkToHead(), which really needs to
   * specify access level 'private' to ensure non-pipeline classes cannot
   * accidentally invoke it.
   * Another solution would be to designate SubmitThreadWork() as a virtual
   * method, and then override it in AbstractPipelineHead, however this would
   * enable classes other than AbstractPipelineHead to override it as well,
   * potentially breaking the design.
   */
  friend class AbstractPipelineHead;

public:

  AbstractPipelineNode ();

  virtual ~AbstractPipelineNode ();

  enum class OutputPreservationMode {
    kNone,
    kReference,
    kCopy
  };

  /** A shared_ptr to the node's most recent output data is preserved if
   * a mode other than OutputPreservationMode::kNone is set;
   * kReference simply stores a shared_ptr to the originally output ProcData,
   * which incurs virtually no additional overhead and should therefore be
   * favored over kCopy unless it is known in advance that a downstream
   * pipeline node will overwrite the original output. For example, by default,
   * a FFTW complex-to-real transform utilizes its input buffer while computing
   * its output, wrecking the original input content.
   */
  OutputPreservationMode output_preservation_mode () const {
    return output_preservation_mode_;
  }

  void set_output_preservation_mode (OutputPreservationMode value) {
    output_preservation_mode_ = value;
  }

  /** It is the caller's responsibility to ensure thread safety when invoking
   * this getter, as otherwise it may be possible for the returned ProcData
   * to be modified while in use by the caller.
   * Further, the returned shared_ptr will be empty unless a mode other
   * than OutputPreservationMode::kNone has been configured.
   */
  std::shared_ptr<const ProcData> most_recent_output () const {
    return most_recent_output_;
  }

  const std::string& process_error_str () const {
    return process_error_str_;
  }

  void set_node_debug_enabled (bool value) {
    node_debug_enabled_ = value;
  }

  std::string ConnectOutput (std::shared_ptr<AbstractPipelineNode> output_node);

protected:

  std::shared_ptr<AbstractPipelineNode> output_node () {
    return output_node_;
  }

  bool node_debug_enabled () const {
    return node_debug_enabled_;
  }

  void SubmitThreadWork (ThreadSquad::ThreadWorkUnit work);

  enum class ReturnReason { kTimeExpired, kShuttingDown };

  /** See documentation for AbstractPipelineHead::InterruptibleSleep()
   */
  ReturnReason InterruptibleSleep (double seconds);

private:

  /** Intended to be accessed exclusively by AbstractPipelineHead
   */
  void set_most_recent_output (std::shared_ptr<ProcData> value) {
    most_recent_output_ = value;
  }

  void set_process_error_str (const std::string& value) {
    process_error_str_ = value;
  }

  /** Naming convention IsCamelCase to distinguish it from trivial_setters;
   * this method does more than simply reassign an internal var
   */
  void SetHead (AbstractPipelineHead* head);

  /** Invoked exclusively by this->head_ during destruction to nullify all
   * pointers to the head object
   */
  void InvalidateHead ();

  virtual std::string Process (std::shared_ptr<const ProcData> input,
                               std::shared_ptr<ProcData>* output) = 0;

  /** Each node's reference to the pipeline head is in the form of a raw pointer
   * instead of a managed shared_ptr because:
   * If a node refers to the head through a shared_ptr, and the head in turn
   * refers to the node through a shared_ptr (i.e., its output_node_), then
   * the reference count of both shared pointers will never drop below zero;
   * i.e., the mutual shared_ptr references will cause the underlying node
   * instances to become lost memory when both nodes go out of scope.
   * Consider also that the head node having a shared_ptr to itself would
   * prevent a head node instance from being automatically destructed.
   */
  AbstractPipelineNode* head_;

  /** Node-specific debug setting.
   * See also AbstractPipelineHead::pipeline_debug_enabled_
   */
  bool node_debug_enabled_;

  std::shared_ptr<AbstractPipelineNode> output_node_;

  /** A shared_ptr to the node's most recent output data is preserved if
   * a mode other than OutputPreservationMode::kNone is set
   */
  OutputPreservationMode output_preservation_mode_;

  std::shared_ptr<ProcData> most_recent_output_;

  std::string process_error_str_;
};

}
