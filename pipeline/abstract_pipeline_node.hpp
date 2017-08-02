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

  std::string ConnectOutput (std::shared_ptr<AbstractPipelineNode> output_node);

protected:

  std::shared_ptr<AbstractPipelineNode> output_node () {
    return output_node_;
  }

  void SubmitThreadWork (ThreadSquad::ThreadWorkUnit work);

  enum class ReturnReason { kTimeExpired, kShuttingDown };

  /** See documentation for AbstractPipelineHead::InterruptibleSleep()
   */
  ReturnReason InterruptibleSleep (double seconds);

private:

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

  std::shared_ptr<AbstractPipelineNode> output_node_;
};

}
